"""Byte count for a pointer stitch versus a rounded block read.

Each sample is stored once. A child records parent ids and does not copy
their columns. ``stitch_bytes`` sums the unique columns a trio names.
``rounded_read`` opens every fixed-size block those columns sit in and
counts the whole block, including a short tail.

This is the difference a controller would stream if it followed the
pointers. No CXL pool, SmartNIC, or storage controller was attached.
"""

from __future__ import annotations


class PointerPool:
    """One column per sample, packed in arrival order."""

    def __init__(self, block_bytes: int, sample_bytes: int) -> None:
        if sample_bytes < 1 or block_bytes < sample_bytes:
            raise ValueError("a block must hold at least one sample")
        self.block_bytes = block_bytes
        self.sample_bytes = sample_bytes
        self.slots: list[int] = []
        self.parents: dict[int, tuple[int, int]] = {}

    def add_sample(self, sample_id: int, father: int = -1, mother: int = -1) -> None:
        if sample_id in self._index():
            raise ValueError(f"sample {sample_id} is already stored")
        if father >= 0 or mother >= 0:
            if father < 0 or mother < 0:
                raise ValueError("a child needs both parents")
            index = self._index()
            if father not in index or mother not in index:
                raise ValueError("both parents must already be stored")
            self.parents[sample_id] = (father, mother)
        self.slots.append(sample_id)

    def stored_bytes(self) -> int:
        return len(self.slots) * self.sample_bytes

    def stitch_bytes(self, sample_ids: list[int]) -> int:
        """Unique column bytes. Shared parents are counted once."""

        seen: set[int] = set()
        index = self._index()
        for sample_id in sample_ids:
            if sample_id not in index:
                raise ValueError(f"sample {sample_id} is not stored")
            seen.add(sample_id)
        return len(seen) * self.sample_bytes

    def rounded_read(self, sample_ids: list[int]) -> int:
        """Full blocks touched by those columns."""

        wanted = set(sample_ids)
        opened = 0
        found = set()
        for block in self._packed():
            hit = [sample_id for sample_id in block if sample_id in wanted]
            if hit:
                opened += self.block_bytes
                found.update(hit)
        if found != wanted:
            raise ValueError("query samples are missing from the pool")
        return opened

    def materialized_bytes(self, trios: list[tuple[int, int, int]]) -> int:
        """Bytes if every trio is copied into its own buffer."""

        total = 0
        index = self._index()
        for father, mother, child in trios:
            for sample_id in (father, mother, child):
                if sample_id not in index:
                    raise ValueError(f"sample {sample_id} is not stored")
            total += 3 * self.sample_bytes
        return total

    def _index(self) -> dict[int, int]:
        return {sample_id: position for position, sample_id in enumerate(self.slots)}

    def _packed(self) -> list[list[int]]:
        blocks: list[list[int]] = []
        current: list[int] = []
        used = 0
        for sample_id in self.slots:
            if current and used + self.sample_bytes > self.block_bytes:
                blocks.append(current)
                current = []
                used = 0
            current.append(sample_id)
            used += self.sample_bytes
        if current:
            blocks.append(current)
        return blocks
