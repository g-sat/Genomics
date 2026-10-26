# Trio accumulator

`src/pedpack/trio_accumulator.py`. Synthetic ids only. The planner holds the HMAC key and the parent links. The store receives a Merkle root and one opaque payload per trio leaf. The leaf is the hash of three sorted tokens. A fetch is accepted only when the sibling path recomputes that root.

Thirty-two trios, 120-byte payloads:

| Item | Value |
| --- | ---: |
| Store document | 10072 bytes |
| Membership proof | 192 bytes on every trio |
| Bytes returned for one trio | 120 |
| Relationship words in the document | none |
| HMAC token present in the document | no |

A sibling with one flipped bit is rejected, and so is a leaf built with a different key. The other thirty-one payloads stay unopened.

The store sees the leaf that was opened. It does not see sample ids or parent edges. Dosages are not computed under encryption, and this path does not run in an enclave. The planner that holds the key can recompute every token. This is not a HIPAA or GDPR control.
