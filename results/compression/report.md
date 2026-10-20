# Zstandard on nuclear-family tiles

Simulated haplotypes and dosages, same generator and seeds as the studies.
Codec is Zstandard level 3. Each tile is 800 variants,
variant-major, and the control is an equal-sized draw of founders.
The model ratio is `payload(family kinship) / payload(kinship 0)` with
`compress_base=0.45` and `kinship_alpha=0.55`.
This is not PBWT. A generic codec can miss a haplotype transform that
PBWT would catch, and the other way around.

| Regime | People | Groups | Mean kinship | Hap family/raw | Hap founder/raw | Hap family/founder | Dosage family/founder | Model family/unrelated |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| locked | 76, 58, 67 | 55 | 0.213 | 0.165 | 0.182 | 0.907 | 0.892 | 0.883 |
| scale | 1207, 1175, 1176 | 1057 | 0.203 | 0.167 | 0.182 | 0.918 | 0.903 | 0.888 |
