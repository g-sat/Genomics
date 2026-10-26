"""Token index for a pedigree the block store does not store.

The planner keeps parent pointers. The store receives an HMAC token per
sample and the block that token lives in. A hereditary check is a list of
tokens. This is aliasing with a secret key. It is not a zero-knowledge
proof and it is not homomorphic encryption: anyone who has the key and a
sample id can recompute that token.
"""

from __future__ import annotations

import hashlib
import hmac
import json


def sample_token(key: bytes, sample_id: int) -> str:
    if not key:
        raise ValueError("the planner key must be non-empty")
    digest = hmac.new(key, str(sample_id).encode("ascii"), hashlib.sha256)
    return digest.hexdigest()[:32]


class TokenIndex:
    """Block map published to storage. Keys are tokens, never sample ids."""

    def __init__(self, key: bytes) -> None:
        self._key = key
        self.block_of: dict[str, int] = {}

    def publish(self, block_id: int, sample_ids: list[int]) -> None:
        for sample_id in sample_ids:
            self.block_of[sample_token(self._key, sample_id)] = block_id

    def plan(self, sample_ids: list[int]) -> list[str]:
        return [sample_token(self._key, sample_id) for sample_id in sample_ids]

    def blocks_for(self, tokens: list[str]) -> list[int]:
        return [self.block_of[token] for token in tokens]

    def storage_document(self) -> str:
        """JSON the store is allowed to see."""

        return json.dumps({"blocks": self.block_of}, sort_keys=True)


_PEDIGREE_WORDS = ("father", "mother", "parent", "kinship", "trio", "pedigree")


def document_names_a_relationship(document: str) -> bool:
    lowered = document.lower()
    return any(word in lowered for word in _PEDIGREE_WORDS)
