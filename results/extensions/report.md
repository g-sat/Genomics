# Three layout extensions

Synthetic samples only. No biobank, no clinical notes, and no encrypted-query proof.

## Append-only pedigree

`src/pedpack/incremental.py`. Block 200 bytes, sample 40 bytes, six samples, two trios.

| Event | Payload rewritten | Blocks appended | Full rebuild |
| --- | ---: | ---: | ---: |
| After the six inserts | 0 | 200 | 400 |
| After sealing trio (0, 1, 2) | 120 | 600 | 400 |

The seal copies three samples. The other three stay in the block they already occupied. Parent bytes are not rewritten when the child is appended.

## Token index

`src/pedpack/private_index.py`. The store's JSON is an HMAC token mapped to a block id. The document does not contain father, mother, parent, kinship, trio, or pedigree. A check is three tokens. The same id under a different key is a different token.

This hides the edge list from the store. It is not a zero-knowledge proof and it is not homomorphic encryption. The planner that holds the key can recompute every token.

## Pedigree cut versus phenotype cut

`src/pedpack/multimodal.py`. Nine synthetic ids. Pedigree groups are the three trios. Phenotype groups take one person from each trio. A hereditary check of one trio reads 300 bytes on the pedigree cut and 900 on the phenotype cut. The swapped check reads 900 and 300. The pedigree cut is cheaper when more than half of the queries are hereditary, and the phenotype cut is cheaper below that. At one half the two costs tie.
