# Pointer stitch

`src/pedpack/pointer_resolve.py`. Each sample is stored once. A child records parent ids and does not copy their columns. The stitch sums unique column bytes. The block path charges every 64 KB block those columns touch, including a short tail. Column width is 200 bytes, the same window width as the public-callset check. No CXL pool, SmartNIC, or storage controller was attached.

| Read | Rounded blocks | Pointer stitch |
| --- | ---: | ---: |
| One trio packed in the first block | 65536 | 600 |
| Two trios that share both parents | 65536 | 800 |
| Same trio with one person in each of three blocks | 196608 | 600 |

The two shared parents occupy 800 stored bytes. Copying each trio into its own buffer writes 1200 bytes. The stitch of both children still returns 800, and it does not allocate that extra copy.

On the small blocks in `tests/test_pointer_resolve.py` the same gap appears at a finer scale: a 3-column read rounds to 512 bytes against a 300-byte stitch, and three columns in separate 200-byte blocks round to 600 against 300.
