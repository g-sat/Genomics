"""Merkle accumulator for trio blocks the store cannot name.

The planner keeps the sample ids and the HMAC key. It publishes a Merkle
root and one opaque payload per trio commitment. A fetch presents the leaf
and the sibling path. The store returns that payload only when the path
recomputes the root.

The leaf is the hash of three sorted HMAC tokens, so the published document
has no sample id and no parent edge. A valid path still reveals the leaf
that was opened. This is an accumulator check, not a zero-knowledge proof,
not homomorphic encryption, and not an enclave. Anyone who holds the key
can recompute a token.
"""

from __future__ import annotations

import hashlib
import hmac
import json

from pedpack.private_index import sample_token

_EMPTY = b"\x00" * 32


def _parent(left: bytes, right: bytes) -> bytes:
    return hashlib.sha256(b"\x01" + left + right).digest()


def merkle_root_and_proofs(
    leaves: list[bytes],
) -> tuple[bytes, list[list[tuple[bytes, bool]]]]:
    """Sibling path for every leaf. ``True`` means the sibling is on the left."""

    if not leaves:
        raise ValueError("the accumulator needs at least one trio")
    proofs: list[list[tuple[bytes, bool]]] = [[] for _ in leaves]
    level = list(leaves)
    owners: list[list[int]] = [[index] for index in range(len(leaves))]
    while len(level) > 1:
        nxt: list[bytes] = []
        nxt_owners: list[list[int]] = []
        index = 0
        while index < len(level):
            if index + 1 == len(level):
                parent = _parent(level[index], _EMPTY)
                for owner in owners[index]:
                    proofs[owner].append((_EMPTY, False))
                nxt.append(parent)
                nxt_owners.append(owners[index])
                index += 1
                continue
            parent = _parent(level[index], level[index + 1])
            for owner in owners[index]:
                proofs[owner].append((level[index + 1], False))
            for owner in owners[index + 1]:
                proofs[owner].append((level[index], True))
            nxt.append(parent)
            nxt_owners.append(owners[index] + owners[index + 1])
            index += 2
        level = nxt
        owners = nxt_owners
    return level[0], proofs


def verify_membership(
    root: bytes,
    leaf: bytes,
    proof: list[tuple[bytes, bool]],
) -> bool:
    digest = leaf
    for sibling, sibling_is_left in proof:
        if sibling_is_left:
            digest = _parent(sibling, digest)
        else:
            digest = _parent(digest, sibling)
    return hmac.compare_digest(digest, root)


def open_block(
    root: bytes,
    blocks: dict[str, str],
    leaf: bytes,
    proof: list[tuple[bytes, bool]],
) -> bytes:
    """Storage-node fetch. The arguments are hashes, not sample ids."""

    if not verify_membership(root, leaf, proof):
        raise PermissionError("accumulator rejected the proof")
    payload = blocks.get(leaf.hex())
    if payload is None:
        raise PermissionError("accumulator rejected the proof")
    return bytes.fromhex(payload)


class TrioAccumulator:
    """Client-side planner. ``storage_document`` is the only store input."""

    def __init__(self, key: bytes) -> None:
        if not key:
            raise ValueError("the planner key must be non-empty")
        self._key = key
        self.root = b""
        self._blocks: dict[str, str] = {}
        self._proofs: dict[bytes, list[tuple[bytes, bool]]] = {}
        self._leaves: dict[tuple[int, int, int], bytes] = {}

    def leaf_for(self, members: tuple[int, int, int]) -> bytes:
        tokens = sorted(sample_token(self._key, sample_id) for sample_id in members)
        raw = "\n".join(tokens).encode("ascii")
        return hashlib.sha256(b"commitment\n" + raw).digest()

    def publish(
        self,
        trios: list[tuple[int, int, int]],
        payloads: list[bytes],
    ) -> None:
        if len(trios) != len(payloads):
            raise ValueError("each trio needs one payload")
        leaves: list[bytes] = []
        blocks: dict[str, str] = {}
        known: dict[tuple[int, int, int], bytes] = {}
        for trio, payload in zip(trios, payloads, strict=True):
            leaf = self.leaf_for(trio)
            if leaf.hex() in blocks:
                raise ValueError("duplicate trio commitment")
            leaves.append(leaf)
            blocks[leaf.hex()] = payload.hex()
            known[trio] = leaf
        root, proofs = merkle_root_and_proofs(leaves)
        self.root = root
        self._blocks = blocks
        self._leaves = known
        self._proofs = {leaf: proof for leaf, proof in zip(leaves, proofs, strict=True)}

    def storage_document(self) -> str:
        body = {"root": self.root.hex(), "blocks": self._blocks}
        return json.dumps(body, sort_keys=True)

    def proof_for(
        self, members: tuple[int, int, int]
    ) -> tuple[bytes, list[tuple[bytes, bool]]]:
        leaf = self._leaves.get(members)
        if leaf is None:
            raise PermissionError("accumulator rejected the proof")
        return leaf, self._proofs[leaf]

    def fetch(self, members: tuple[int, int, int]) -> bytes:
        leaf, proof = self.proof_for(members)
        return open_block(self.root, self._blocks, leaf, proof)

    def proof_bytes(self, members: tuple[int, int, int]) -> int:
        leaf, proof = self.proof_for(members)
        return len(leaf) + sum(len(sibling) for sibling, _side in proof)
