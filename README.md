# vcf-qc-toolkit

A dependency-free Python tool that summarises the quality of a VCF file in one
pass: variant types, transition/transversion ratio, FILTER status, QUAL,
allele-frequency spectrum, indel length distribution, and per-sample genotype
counts. Reads plain or gzipped VCF using only the standard library.

These are the first sanity checks usually run on a new call set, for example
before or after filtering with `bcftools`.

## Install and run

```bash
pip install -e .
python -m vcfqc examples/example.vcf          # text report
python -m vcfqc calls.vcf.gz --json > qc.json # machine-readable
```

Example output (`examples/example.vcf`):

```text
Records:            6
Variant types:      DEL=1, INS=1, SNP=5
Ti/Tv (SNPs):       4.000  (4 ts / 1 tv)
FILTER:             FAIL=1, PASS=5
QUAL mean/median:   35.00 / 35.00
Per chromosome:     chr1=3, chr2=3
Indel lengths:      -1:1, +2:1
AF histogram:       1 1 0 1 0 1 0 0 0 1

Sample              hom_ref   het  hom_alt  missing  het/hom_alt
S1                       1     4        1        0         4.00
S2                       2     2        1        1         2.00
```

## Python API

```python
from vcfqc import parse_vcf, summarize

with open("calls.vcf") as fh:
    samples, variants = parse_vcf(fh)
    summary = summarize(samples, variants)
print(summary["ti_tv"])
```

## How the metrics are defined

- **Counting unit:** variant types and Ti/Tv are counted **per alternate
  allele**, so a multiallelic record contributes one count per ALT. `n_records`
  counts VCF lines.
- **Variant types:** `SNP`, `MNP` (equal length > 1), `INS`, `DEL`,
  `SYMBOLIC` (`<DEL>`, `*`, breakend notation) and `OTHER` (for example `N`).
  Indel lengths are `len(ALT) - len(REF)`.
- **Ti/Tv:** transitions (A<->G, C<->T) over transversions, SNPs only. For
  whole-genome human data this is typically around 2.0-2.1 and around 3 for
  exomes; values far from the expected range can indicate false positives.
  Small or targeted call sets can legitimately differ.
- **FILTER:** `PASS` and `.` count as pass; anything else counts as fail.
- **AF histogram:** the `AF` INFO field (first value for multiallelic
  records) in 10 equal bins from 0 to 1; records without a valid `AF` are
  skipped.
- **Genotypes:** `het` means the called alleles differ (so `1/2` is het);
  any missing allele makes the call `missing`. `het_hom_alt_ratio` is `None`
  when a sample has no hom-alt calls. Phasing (`|`) is ignored.

## Limitations

- Reads the whole file in a single streaming pass but keeps QUAL values in
  memory for the median; this is fine for typical files, though very large
  multi-sample VCFs will be slow in pure Python.
- Only `GT` is read from the sample columns; depth, genotype quality and other
  FORMAT fields are not summarised.
- No reference-aware checks (for example REF matching the genome).

## Tests

```bash
pip install pytest
pytest
```
