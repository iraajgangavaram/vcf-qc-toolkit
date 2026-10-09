import pytest

from vcfqc import classify, genotype_class, is_transition, parse_vcf, summarize

HEADER = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"


def _vcf(*records, samples=("S1", "S2")):
    cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]
    lines = ["\t".join(cols)] if samples else ["\t".join(cols[:8])]
    lines += ["\t".join(r) for r in records]
    return (HEADER + "\n".join(lines) + "\n").splitlines(keepends=True)


# Every expected value in test_summary_matches_hand_computation was worked out
# by hand from these six records.
RECORDS = [
    ("chr1", "100", ".", "A", "G", "50", "PASS", "AF=0.5", "GT", "0/1", "1/1"),
    ("chr1", "200", ".", "C", "T", "40", "PASS", "AF=0.1", "GT", "0/0", "0/1"),
    ("chr1", "300", ".", "A", "C", "30", "q10", "AF=0.9", "GT", "0|1", "./."),
    ("chr2", "50", ".", "G", "A", "60", "PASS", "AF=0.02", "GT", "1/1", "0/1"),
    ("chr2", "60", ".", "AT", "A", "20", "PASS", "AF=0.3", "GT", "0/1", "0/0"),
    ("chr2", "70", ".", "A", "ATT,G", "10", "PASS", ".", "GT", "1/2", "0/0"),
]


def test_summary_matches_hand_computation():
    samples, variants = parse_vcf(_vcf(*RECORDS))
    s = summarize(samples, variants)

    assert s["n_records"] == 6
    # alt alleles: 4 SNPs + 1 SNP from the multiallelic record, 1 insertion, 1 deletion
    assert s["variant_types"] == {"SNP": 5, "INS": 1, "DEL": 1}
    # transitions: A>G, C>T, G>A, A>G; transversion: A>C
    assert (s["n_transitions"], s["n_transversions"]) == (4, 1)
    assert s["ti_tv"] == 4.0
    assert s["per_chromosome"] == {"chr1": 3, "chr2": 3}
    assert s["filter"] == {"PASS": 5, "FAIL": 1}
    assert s["qual_mean"] == 35.0 and s["qual_median"] == 35.0
    assert s["indel_length_distribution"] == {-1: 1, 2: 1}
    # AFs 0.5, 0.1, 0.9, 0.02, 0.3 (record 6 has none) -> bins 5, 1, 9, 0, 3
    assert s["af_histogram"] == [1, 1, 0, 1, 0, 1, 0, 0, 0, 1]
    # S1: het = 0/1, 0|1, 0/1, 1/2; hom_ref = 0/0; hom_alt = 1/1
    assert s["samples"]["S1"] == {
        "hom_ref": 1, "het": 4, "hom_alt": 1, "missing": 0, "het_hom_alt_ratio": 4.0,
    }
    # S2: hom_alt 1/1; het 0/1 x2; missing ./.; hom_ref 0/0 x2
    assert s["samples"]["S2"] == {
        "hom_ref": 2, "het": 2, "hom_alt": 1, "missing": 1, "het_hom_alt_ratio": 2.0,
    }


@pytest.mark.parametrize(
    "ref,alt,expected",
    [
        ("A", "G", "SNP"),
        ("AT", "GC", "MNP"),
        ("A", "ATT", "INS"),
        ("AT", "A", "DEL"),
        ("A", "<DEL>", "SYMBOLIC"),
        ("A", "*", "SYMBOLIC"),
        ("A", "G]chr2:5]", "SYMBOLIC"),
        ("A", "N", "OTHER"),
    ],
)
def test_classify(ref, alt, expected):
    assert classify(ref, alt) == expected


def test_transition_definition():
    assert is_transition("A", "G") and is_transition("G", "A")
    assert is_transition("C", "T") and is_transition("T", "C")
    assert not is_transition("A", "C") and not is_transition("G", "T")


@pytest.mark.parametrize(
    "gt,expected",
    [
        ("0/0", "hom_ref"), ("0|0", "hom_ref"), ("1/1", "hom_alt"), ("2/2", "hom_alt"),
        ("0/1", "het"), ("1|0", "het"), ("1/2", "het"),
        ("./.", "missing"), ("./1", "missing"), (".", "missing"),
    ],
)
def test_genotype_class(gt, expected):
    assert genotype_class(gt) == expected


def test_vcf_without_samples_or_quality_or_af():
    rec = ("chr1", "5", ".", "A", "G", ".", ".", ".")
    samples, variants = parse_vcf(_vcf(rec, samples=()))
    s = summarize(samples, variants)
    assert samples == [] and s["samples"] == {}
    assert s["qual_mean"] is None and s["filter"] == {"PASS": 1}
    assert sum(s["af_histogram"]) == 0


def test_af_edge_values_and_multiallelic_af_uses_first_value():
    recs = [
        ("chr1", "1", ".", "A", "G", "9", "PASS", "AF=1.0", "GT", "1/1", "1/1"),
        ("chr1", "2", ".", "A", "G,C", "9", "PASS", "AF=0.25,0.05", "GT", "0/1", "0/1"),
        ("chr1", "3", ".", "A", "G", "9", "PASS", "AF=oops", "GT", "0/1", "0/1"),
    ]
    samples, variants = parse_vcf(_vcf(*recs))
    hist = summarize(samples, variants)["af_histogram"]
    assert hist[9] == 1 and hist[2] == 1 and sum(hist) == 2


def test_het_ratio_is_none_when_no_hom_alt():
    rec = ("chr1", "1", ".", "A", "G", "9", "PASS", ".", "GT", "0/1", "0/0")
    samples, variants = parse_vcf(_vcf(rec))
    assert summarize(samples, variants)["samples"]["S1"]["het_hom_alt_ratio"] is None


def test_malformed_input_raises():
    with pytest.raises(ValueError):
        parse_vcf(["##fileformat=VCFv4.2\n"])
    with pytest.raises(ValueError):
        parse_vcf(["chr1\t1\t.\tA\tG\n"])
    _, variants = parse_vcf(_vcf(("chr1", "x", ".", "A", "G", "9", "PASS", ".", "GT", "0/1", "0/1")))
    with pytest.raises(ValueError):
        list(variants)
    _, variants = parse_vcf(_vcf(("chr1", "1", ".", "A")))
    with pytest.raises(ValueError):
        list(variants)
