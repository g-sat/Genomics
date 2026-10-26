"""Append-only pedigree layout with family compaction.

New samples land in a memory tail. A full tail is sealed into an immutable
level-0 run. A hereditary check follows the alias of each person and does
not scan every run. When a nuclear family is spread over more than one
block, compaction appends a level-1 run that holds only that family and
swings the aliases. Older runs stay on the log.

This is the layout rule, not a disk engine. It does not read a biobank.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class _Run:
    level: int
    blocks: list[tuple[int, ...]] = field(default_factory=list)


class PedigreeLSM:
    """Two-level log: arrival-order runs, then one run per compacted family."""

    def __init__(self, block_bytes: int, sample_bytes: int) -> None:
        if block_bytes < sample_bytes or sample_bytes < 1:
            raise ValueError("a block must hold at least one sample")
        self.block_bytes = block_bytes
        self.sample_bytes = sample_bytes
        self._tail: list[int] = []
        self._tail_used = 0
        self.runs: list[_Run] = []
        self.alias: dict[int, tuple[int, int]] = {}
        self.parents: dict[int, tuple[int, int]] = {}
        self.bytes_written = 0
        self.compaction_payload = 0
        self.largest_compaction = 0

    def add_sample(self, sample_id: int, father: int = -1, mother: int = -1) -> None:
        """Append one genome. Parents already stored are not copied."""

        if sample_id in self.alias or sample_id in self._tail:
            raise ValueError(f"sample {sample_id} is already stored")
        if self._tail and self._tail_used + self.sample_bytes > self.block_bytes:
            self.flush()
        self._tail.append(sample_id)
        self._tail_used += self.sample_bytes
        if father >= 0 or mother >= 0:
            if father < 0 or mother < 0:
                raise ValueError("a child needs both parents")
            if not self._stored(father) or not self._stored(mother):
                raise ValueError("both parents must already be stored")
            self.parents[sample_id] = (father, mother)

    def flush(self) -> None:
        """Seal the tail as one level-0 run. Live samples are not rewritten."""

        if not self._tail:
            return
        run = _Run(level=0)
        self.runs.append(run)
        run_index = len(self.runs) - 1
        buf: list[int] = []
        used = 0
        for sample_id in self._tail:
            if buf and used + self.sample_bytes > self.block_bytes:
                self._seal_block(run, run_index, buf)
                buf = []
                used = 0
            buf.append(sample_id)
            used += self.sample_bytes
        if buf:
            self._seal_block(run, run_index, buf)
        self._tail = []
        self._tail_used = 0

    def family_blocks(self, child: int) -> set[tuple[int, int]]:
        """Locations a trio check opens. The open tail is run -1."""

        father, mother = self.parents[child]
        return {self._locate(father), self._locate(mother), self._locate(child)}

    def compact_family(self, child: int) -> int:
        """Append a level-1 run for one nuclear family. Return blocks written.

        The payload counter grows by the family only. Every other alias stays.
        """

        self.flush()
        father, mother = self.parents[child]
        members = [father, mother]
        for other, pair in self.parents.items():
            if pair == (father, mother):
                members.append(other)
        payload = len(members) * self.sample_bytes
        run = _Run(level=1)
        self.runs.append(run)
        run_index = len(self.runs) - 1
        written = 0
        buf: list[int] = []
        used = 0
        for sample_id in members:
            if buf and used + self.sample_bytes > self.block_bytes:
                self._seal_block(run, run_index, buf)
                written += 1
                buf = []
                used = 0
            buf.append(sample_id)
            used += self.sample_bytes
        if buf:
            self._seal_block(run, run_index, buf)
            written += 1
        self.compaction_payload += payload
        self.largest_compaction = max(self.largest_compaction, payload)
        return written

    def compact_scattered(self) -> int:
        """Compact each family whose live copies occupy more than one block."""

        families = 0
        for child in list(self.parents):
            if len(self.family_blocks(child)) > 1:
                self.compact_family(child)
                families += 1
        return families

    def full_rebuild_bytes(self) -> int:
        """One packed rewrite of every person stored so far."""

        n = len(self.alias) + len(self._tail)
        if n == 0:
            return 0
        per_block = self.block_bytes // self.sample_bytes
        blocks = (n + per_block - 1) // per_block
        return blocks * self.block_bytes

    def level_counts(self) -> dict[int, int]:
        counts: dict[int, int] = {}
        for run in self.runs:
            counts[run.level] = counts.get(run.level, 0) + 1
        return counts

    def _stored(self, sample_id: int) -> bool:
        return sample_id in self.alias or sample_id in self._tail

    def _locate(self, sample_id: int) -> tuple[int, int]:
        if sample_id in self._tail:
            return (-1, 0)
        return self.alias[sample_id]

    def _seal_block(self, run: _Run, run_index: int, sample_ids: list[int]) -> None:
        block_index = len(run.blocks)
        run.blocks.append(tuple(sample_ids))
        for sample_id in sample_ids:
            self.alias[sample_id] = (run_index, block_index)
        self.bytes_written += self.block_bytes
