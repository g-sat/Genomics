"""Pointer stitch, deletion aliases, and Parquet row groups on real chr22 dosages.

Reads the high-coverage chromosome 22 panel. Dosages are counted GT alleles.
The packed file is sample-major, one column per sample, padded to 64 KB blocks.
The stitch is the three columns of one real trio with no padding.

    python experiments/real_stitch.py
"""

from __future__ import annotations

import ctypes
import json
import statistics
import sys
import time
from ctypes import wintypes
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import real_trios  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "real_stitch"
BLOCK = 65536
ALIGN = 4096
TRIALS = 7


def _pack(matrix: np.ndarray, block: int) -> tuple[bytes, list[int]]:
    """Sample-major bytes and the block index of each sample. Sites are rows."""

    n_sites, n_samples = matrix.shape
    if n_sites > block:
        raise ValueError("a sample column must fit in one block")
    per = block // n_sites
    chunks: list[bytes] = []
    block_of: list[int] = []
    for start in range(0, n_samples, per):
        piece = bytearray(block)
        used = 0
        for sample in range(start, min(start + per, n_samples)):
            column = matrix[:, sample].tobytes()
            piece[used : used + n_sites] = column
            used += n_sites
            block_of.append(start // per)
        chunks.append(bytes(piece))
    return b"".join(chunks), block_of


def _write(path: Path, payload: bytes) -> None:
    path.write_bytes(payload)


def _open_direct(path: Path) -> int:
    kernel = ctypes.windll.kernel32
    kernel.CreateFileW.argtypes = [
        wintypes.LPCWSTR,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.LPVOID,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.HANDLE,
    ]
    kernel.CreateFileW.restype = wintypes.HANDLE
    handle = kernel.CreateFileW(
        str(path),
        0x80000000,
        0x00000001,
        None,
        3,
        0x20000000,
        None,
    )
    if handle == wintypes.HANDLE(-1).value:
        raise ctypes.WinError()
    return handle


def _read_direct(handle: int, offset: int, size: int) -> int:
    kernel = ctypes.windll.kernel32
    kernel.SetFilePointerEx.argtypes = [
        wintypes.HANDLE,
        ctypes.c_longlong,
        ctypes.POINTER(ctypes.c_longlong),
        wintypes.DWORD,
    ]
    kernel.SetFilePointerEx.restype = wintypes.BOOL
    kernel.ReadFile.argtypes = [
        wintypes.HANDLE,
        wintypes.LPVOID,
        wintypes.DWORD,
        ctypes.POINTER(wintypes.DWORD),
        wintypes.LPVOID,
    ]
    kernel.ReadFile.restype = wintypes.BOOL
    moved = ctypes.c_longlong()
    if not kernel.SetFilePointerEx(handle, offset, ctypes.byref(moved), 0):
        raise ctypes.WinError()
    buffer = ctypes.create_string_buffer(size)
    got = wintypes.DWORD()
    if not kernel.ReadFile(handle, buffer, size, ctypes.byref(got), None):
        raise ctypes.WinError()
    return int(got.value)


def _close(handle: int) -> None:
    ctypes.windll.kernel32.CloseHandle(handle)


def _median_read(path: Path, spans: list[tuple[int, int]]) -> float:
    handle = _open_direct(path)
    try:
        for offset, size in spans:
            _read_direct(handle, offset, size)
        times = []
        for _ in range(TRIALS):
            started = time.perf_counter()
            for offset, size in spans:
                _read_direct(handle, offset, size)
            times.append(time.perf_counter() - started)
    finally:
        _close(handle)
    return statistics.median(times)


def _membership(groups: list[list[int]]) -> dict[int, list[int]]:
    where: dict[int, list[int]] = {}
    for group_index, group in enumerate(groups):
        for sample in group:
            where.setdefault(sample, []).append(group_index)
    return where


def _parquet(matrix: np.ndarray, groups: list[list[int]], trio: list[int], out: Path) -> dict:
    import pyarrow as pa
    import pyarrow.parquet as pq

    n_sites, n_samples = matrix.shape
    dosage = pa.array([matrix[:, sample].tobytes() for sample in range(n_samples)], type=pa.binary())
    table = pa.table({"dosage": dosage})
    whole = out / "cohort.parquet"
    pq.write_table(table, whole, compression="NONE")
    grouped = out / "families.parquet"
    writer = pq.ParquetWriter(grouped, table.schema, compression="NONE")
    trio_group = None
    for index, group in enumerate(groups):
        if set(trio).issubset(group):
            trio_group = index
        writer.write_table(table.take(group))
    writer.close()
    if trio_group is None:
        raise RuntimeError("the trio is not inside one family row group")
    whole_meta = pq.ParquetFile(whole).metadata
    grouped_meta = pq.ParquetFile(grouped).metadata
    return {
        "whole_file_bytes": whole.stat().st_size,
        "family_file_bytes": grouped.stat().st_size,
        "whole_row_groups": whole_meta.num_row_groups,
        "family_row_groups": grouped_meta.num_row_groups,
        "trio_row_group": trio_group,
        "trio_row_group_bytes": grouped_meta.row_group(trio_group).total_byte_size,
        "samples_in_trio_row_group": len(groups[trio_group]),
    }


def _one_window(matrix: np.ndarray, trio: list[int], groups: list[list[int]]) -> dict:
    payload, block_of = _pack(matrix, BLOCK)
    opened = sorted({block_of[sample] for sample in trio})
    stitch = b"".join(matrix[:, sample].tobytes() for sample in trio)
    cohort_path = OUT / f"cohort_{matrix.shape[0]}.bin"
    stitch_path = OUT / f"stitch_{matrix.shape[0]}.bin"
    _write(cohort_path, payload)
    padded = stitch + bytes((-len(stitch)) % ALIGN)
    _write(stitch_path, padded)
    spans = [(block * BLOCK, BLOCK) for block in opened]
    cohort_seconds = _median_read(cohort_path, spans)
    whole_seconds = _median_read(cohort_path, [(0, len(payload))])
    stitch_seconds = _median_read(stitch_path, [(0, len(padded))])
    where = _membership(groups)
    busiest = max(where, key=lambda sample: len(where[sample]))
    withdrawn = trio[0]
    return {
        "sites": int(matrix.shape[0]),
        "sample_bytes": int(matrix.shape[0]),
        "samples": int(matrix.shape[1]),
        "cohort_file_bytes": len(payload),
        "blocks_in_file": len(payload) // BLOCK,
        "stitch_bytes": len(stitch),
        "stitch_aligned_bytes": len(padded),
        "blocks_opened": len(opened),
        "rounded_bytes": len(opened) * BLOCK,
        "unbuffered_trio_seconds": cohort_seconds,
        "unbuffered_cohort_seconds": whole_seconds,
        "unbuffered_stitch_seconds": stitch_seconds,
        "io": "FILE_FLAG_NO_BUFFERING",
        "trials": TRIALS,
        "withdraw_one_block_bytes": BLOCK,
        "withdraw_full_file_bytes": len(payload),
        "withdraw_sample_groups": len(where[withdrawn]),
        "busiest_sample_groups": len(where[busiest]),
        "people_in_two_groups": sum(1 for groups_of in where.values() if len(groups_of) > 1),
    }


def main() -> None:
    if not real_trios.VCF.exists():
        raise SystemExit(f"missing {real_trios.VCF}")
    OUT.mkdir(parents=True, exist_ok=True)
    print("reading dosages", flush=True)
    sample_ids, matrix, scanned = real_trios.load_biallelic(real_trios.N_VARIANTS)
    trios = real_trios.trios_in_panel(set(sample_ids), real_trios.PED)
    groups = real_trios.family_groups(sample_ids, trios)
    multi = [group for group in groups if len(group) >= 3]
    trio = multi[0]
    window = matrix[: real_trios.WINDOW]
    print("packing 200-site window", flush=True)
    small = _one_window(window, trio, groups)
    print("packing 3000-site matrix", flush=True)
    full = _one_window(matrix, trio, groups)
    print("writing parquet row groups", flush=True)
    parquet = _parquet(matrix, groups, trio, OUT)
    record = {
        "samples": int(matrix.shape[1]),
        "biallelic_sites": int(matrix.shape[0]),
        "variant_rows_scanned": scanned,
        "trios": len(trios),
        "first_trio_size": len(trio),
        "block": BLOCK,
        "window_200": small,
        "sites_3000": full,
        "parquet": parquet,
    }
    (OUT / "measurements.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record, indent=2), flush=True)


if __name__ == "__main__":
    main()
