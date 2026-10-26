"""Append-only pedigree layout.

Genotype bytes already sealed in a block are not rewritten when a child
arrives. The child is packed into the open tail. A trio query follows
aliases to the blocks that hold the father, the mother, and the child.
Sealing a family appends one new copy of that trio and leaves the old
blocks in place. The rewrite counter counts only those copied bytes.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class AppendOnlyLayout:
    block_bytes: int
    sample_bytes: int
    blocks: list[tuple[int, ...]] = field(default_factory=list)
    tail: list[int] = field(default_factory=list)
    tail_used: int = 0
    alias: dict[int, int] = field(default_factory=dict)
    parents: dict[int, tuple[int, int]] = field(default_factory=dict)
    appended_bytes: int = 0
    rewritten_bytes: int = 0

    def __post_init__(self) -> None:
        if self.block_bytes < self.sample_bytes:
            raise ValueError("a block must hold at least one sample")
        if self.sample_bytes < 1:
            raise ValueError("sample payload must be positive")

    def add_sample(self, sample_id: int, father: int = -1, mother: int = -1) -> None:
        """Append one sample. Parents, when given, are not copied."""

        if sample_id in self.alias:
            raise ValueError(f"sample {sample_id} is already stored")
        if self.tail and self.tail_used + self.sample_bytes > self.block_bytes:
            self._seal_tail()
        self.tail.append(sample_id)
        self.tail_used += self.sample_bytes
        self.alias[sample_id] = -1
        if father >= 0 or mother >= 0:
            if father < 0 or mother < 0:
                raise ValueError("a child needs both parents")
            if father not in self.alias or mother not in self.alias:
                raise ValueError("both parents must already be stored")
            self.parents[sample_id] = (father, mother)

    def trio_blocks(self, child: int) -> set[int]:
        """Block ids a hereditary check opens. The tail is block -1."""

        father, mother = self.parents[child]
        return {self.alias[child], self.alias[father], self.alias[mother]}

    def seal_family(self, child: int) -> int:
        """Append a fresh copy of one trio. Return the number of new blocks.

        Old genotype blocks stay on the log. ``rewritten_bytes`` grows by
        the trio payload only.
        """

        father, mother = self.parents[child]
        self._seal_tail()
        copied = 0
        new_blocks = 0
        buf: list[int] = []
        used = 0
        for sample_id in (father, mother, child):
            if buf and used + self.sample_bytes > self.block_bytes:
                self._write(buf)
                new_blocks += 1
                buf = []
                used = 0
            buf.append(sample_id)
            used += self.sample_bytes
            copied += self.sample_bytes
        if buf:
            self._write(buf)
            new_blocks += 1
        self.rewritten_bytes += copied
        return new_blocks

    def full_rebuild_bytes(self) -> int:
        """Bytes to write every live sample again, rounded up to blocks."""

        if not self.alias:
            return 0
        per_block = max(1, self.block_bytes // self.sample_bytes)
        blocks = (len(self.alias) + per_block - 1) // per_block
        return blocks * self.block_bytes

    def _seal_tail(self) -> None:
        if not self.tail:
            return
        self._write(self.tail)
        self.tail = []
        self.tail_used = 0

    def _write(self, sample_ids: list[int]) -> None:
        block_id = len(self.blocks)
        self.blocks.append(tuple(sample_ids))
        for sample_id in sample_ids:
            self.alias[sample_id] = block_id
        self.appended_bytes += self.block_bytes
