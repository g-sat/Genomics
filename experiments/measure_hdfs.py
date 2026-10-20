"""Read materialized layouts through WebHDFS.

This does not invent measurements. If the NameNode is unreachable, it exits
without writing numbers.

    python experiments/measure_hdfs.py --namenode http://127.0.0.1:9870 --root /pedpack

What it records, from the NameNode and the client, for each layout:
bytes returned to the client, wall time, file count, primary length,
space consumed (replication included), and blocks touched according to
GETFILEBLOCKLOCATIONS. A client that is not itself a DataNode sees every
block as remote; that fraction is reported as such and is not the simulator's
task-placement remote fraction.
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


def _with_user(url: str) -> str:
    if "user.name=" in url:
        return url
    joiner = "&" if "?" in url else "?"
    return f"{url}{joiner}user.name=root"


def _request(url: str, data: bytes | None = None, method: str | None = None, timeout: float = 120) -> tuple[int, dict, bytes]:
    request = urllib.request.Request(_with_user(url), data=data, method=method)
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            body = response.read()
            raw_type = response.headers.get("Content-Type", "")
            payload = {}
            if "json" in raw_type or body[:1] in (b"{", b"["):
                try:
                    payload = json.loads(body.decode("utf-8"))
                except json.JSONDecodeError:
                    payload = {}
            return response.status, payload, body
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"WebHDFS {exc.code} for {url}: {detail}") from exc


def _op(namenode: str, path: str, op: str, method: str | None = None, **params: str) -> dict:
    query = {"op": op, **params}
    url = f"{namenode.rstrip('/')}/webhdfs/v1{path}?{urllib.parse.urlencode(query)}"
    _status, payload, _body = _request(url, method=method)
    return payload


def _put(namenode: str, remote: str, blob: bytes) -> None:
    query = urllib.parse.urlencode({"op": "CREATE", "overwrite": "true", "noredirect": "true"})
    url = f"{namenode.rstrip('/')}/webhdfs/v1{remote}?{query}"
    _status, payload, _body = _request(url, method="PUT")
    location = payload.get("Location")
    if not location:
        raise RuntimeError(f"CREATE did not return a Location for {remote}")
    _request(location, data=blob, method="PUT", timeout=600)


def _open(namenode: str, remote: str, offset: int, length: int) -> tuple[int, float]:
    query = urllib.parse.urlencode(
        {"op": "OPEN", "offset": str(offset), "length": str(length), "noredirect": "true"}
    )
    url = f"{namenode.rstrip('/')}/webhdfs/v1{remote}?{query}"
    _status, payload, _body = _request(url)
    location = payload.get("Location")
    if not location:
        raise RuntimeError(f"OPEN did not return a Location for {remote}")
    started = time.perf_counter()
    _status, _payload, body = _request(location, timeout=600)
    elapsed = time.perf_counter() - started
    return len(body), elapsed


def _blocks(namenode: str, remote: str, offset: int, length: int) -> int:
    payload = _op(
        namenode,
        remote,
        "GETFILEBLOCKLOCATIONS",
        offset=str(offset),
        length=str(length),
    )
    locations = payload.get("BlockLocations", {}).get("BlockLocation", [])
    return len(locations)


def measure_layout(namenode: str, root: str, layout_dir: Path) -> dict:
    manifest = json.loads((layout_dir / "manifest.json").read_text(encoding="utf-8"))
    reads = json.loads((layout_dir / "reads.json").read_text(encoding="utf-8"))
    remote_dir = f"{root.rstrip('/')}/{layout_dir.parent.name}/{layout_dir.name}"
    _op(namenode, remote_dir, "MKDIRS", method="PUT")
    for record in manifest["files"]:
        blob = (layout_dir / record["path"]).read_bytes()
        _put(namenode, f"{remote_dir}/{record['path']}", blob)
    summary = _op(namenode, remote_dir, "GETCONTENTSUMMARY")["ContentSummary"]
    by_query: dict[str, dict] = {}
    for item in reads:
        bucket = by_query.setdefault(
            item["query"],
            {"workload": item["workload"], "bytes": 0, "seconds": 0.0, "blocks": 0, "opens": 0},
        )
        path = f"{remote_dir}/{item['file']}"
        nbytes, seconds = _open(namenode, path, int(item["offset"]), int(item["length"]))
        bucket["bytes"] += nbytes
        bucket["seconds"] += seconds
        bucket["blocks"] += _blocks(namenode, path, int(item["offset"]), int(item["length"]))
        bucket["opens"] += 1
    workloads: dict[str, list[dict]] = {}
    for row in by_query.values():
        workloads.setdefault(row["workload"], []).append(row)

    def mean(workload: str, key: str) -> float:
        rows = workloads.get(workload, [])
        if not rows:
            return float("nan")
        return float(sum(row[key] for row in rows) / len(rows))

    return {
        "label": manifest["label"],
        "family": manifest["family"],
        "n_files": int(summary["fileCount"]),
        "primary_bytes": int(summary["length"]),
        "space_consumed_bytes": int(summary["spaceConsumed"]),
        "heldout_family_mean_bytes": mean("heldout_family", "bytes"),
        "heldout_family_mean_seconds": mean("heldout_family", "seconds"),
        "heldout_family_mean_blocks": mean("heldout_family", "blocks"),
        "family_selective_mean_bytes": mean("family_selective", "bytes"),
        "family_selective_mean_seconds": mean("family_selective", "seconds"),
        "family_selective_mean_blocks": mean("family_selective", "blocks"),
        "client_is_not_a_datanode": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Measure materialized layouts through WebHDFS")
    parser.add_argument("--namenode", required=True, help="WebHDFS endpoint, for example http://127.0.0.1:9870")
    parser.add_argument("--root", default="/pedpack")
    parser.add_argument("--transfer", type=Path, default=Path("results/transfer"))
    parser.add_argument("--out", type=Path, default=Path("results/hdfs/measurements.json"))
    parser.add_argument("--scenarios", nargs="*", default=["select_05", "select_15", "mix_50", "overlap"])
    args = parser.parse_args()
    try:
        _op(args.namenode, "/", "GETFILESTATUS")
    except (urllib.error.URLError, RuntimeError, TimeoutError, OSError) as exc:
        raise SystemExit(f"NameNode is not reachable at {args.namenode}: {exc}") from exc
    rows = []
    for scenario in args.scenarios:
        scenario_dir = args.transfer / scenario
        if not scenario_dir.exists():
            continue
        for seed_dir in sorted(scenario_dir.glob("seed_*")):
            for method in ("pedpack", "pedpack_no_pedigree_features", "generic_oracle"):
                layout_dir = seed_dir / method
                if not (layout_dir / "manifest.json").exists():
                    continue
                print(f"measuring {scenario} {seed_dir.name} {method}", flush=True)
                measured = measure_layout(args.namenode, f"{args.root}/{scenario}", layout_dir)
                measured["scenario"] = scenario
                measured["seed_dir"] = seed_dir.name
                measured["method"] = method
                rows.append(measured)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
