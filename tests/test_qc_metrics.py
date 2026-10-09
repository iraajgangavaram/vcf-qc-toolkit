"""Every expected value below was worked out by hand from the five-sample VCF in VCF_TEXT."""

import pytest

from vcfqc import parse_vcf, summarize
from vcfqc.core import robust_z, substitution_class
from vcfqc.simulate import expected_ti_tv, simulate_vcf

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\n"
ROWS = [
    ("chr1", "100", ".", "A", "G", "50", "PASS", ".", "GT:DP:GQ",
     "0/0:10:50", "0/1:20:60", "0/0:30:70", "0/0:40:80", "./.:.:."),
    ("chr1", "200", ".", "C", "A", "50", "PASS", ".", "GT:DP:GQ",
     "1/1:10:40", "0/1:20:50", "0/1:30:60", "0/0:40:70", "0/0:50:80"),
    ("chr2", "50", ".", "G", "A", "50", "PASS", ".", "GT:DP:GQ",
     "0/0:20:30", "0/0:20:30", "0/0:20:30", "0/1:20:30", "0/0:20:30"),
    ("chr2", "60", ".", "T", "C", "50", "PASS", ".", "GT:DP:GQ",
     "1/1:30:90", "1/1:30:90", "1/1:30:90", "1/1:30:90", "0/0:30:90"),
    ("chr3", "10", ".", "A", "C,G", "50", "PASS", ".", "GT:DP:GQ",
     "1/2:.:.", "0/0:.:.", "0/0:.:.", "0/0:.:.", "0/0:.:."),
]
VCF_TEXT = (
    "##fileformat=VCFv4.2\n" + HEADER + "".join("\t".join(r) + "\n" for r in ROWS)
)


def _summary(**kw):
    sites = []
    samples, variants = parse_vcf(VCF_TEXT.splitlines(keepends=True))
    return summarize(samples, variants, site_sink=sites, **kw), sites


def test_substitution_classes_collapse_strands():
    assert substitution_class("A", "G") == "T>C"
    assert substitution_class("G", "A") == "C>T"
    assert substitution_class("T", "G") == "T>G"
    assert substitution_class("C", "A") == "C>A"
    assert substitution_class("A", "T") == "T>A"
    assert substitution_class("G", "C") == "C>G"


def test_substitution_spectrum_and_overall_ti_tv():
    s, _ = _summary()
    # A>G, C>A, G>A, T>C, plus A>C and A>G from the multiallelic record
    assert s["substitution_spectrum"] == {
        "C>A": 1, "C>G": 0, "C>T": 1, "T>A": 0, "T>C": 3, "T>G": 1}
    assert (s["n_transitions"], s["n_transversions"], s["ti_tv"]) == (4, 2, 2.0)


def test_site_level_metrics():
    s, _ = _summary()
    q = s["site_qc"]
    assert q["n_sites_with_genotypes"] == 5
    assert q["mean_site_missingness"] == pytest.approx(0.04)  # only site 1 has 1/5 missing
    assert q["n_sites_missingness_above_threshold"] == 1
    assert (q["singleton_sites"], q["doubleton_sites"]) == (2, 0)
    # observed alt-allele frequencies: .125 (bin 1), .4 (bin 4), .1 (bin 1), .8 (bin 8)
    assert s["af_histogram_observed"] == [0, 2, 0, 0, 1, 0, 0, 0, 1, 0]


def test_ti_tv_by_allele_frequency_bin():
    s, _ = _summary()
    by_bin = s["ti_tv_by_af_bin"]
    assert (by_bin[1]["transitions"], by_bin[1]["transversions"]) == (2, 0)
    assert by_bin[1]["ti_tv"] is None  # no transversions: ratio undefined, not infinity
    assert (by_bin[4]["transitions"], by_bin[4]["transversions"], by_bin[4]["ti_tv"]) == (0, 1, 0.0)
    assert (by_bin[8]["transitions"], by_bin[8]["transversions"]) == (1, 0)


def test_hwe_summary_matches_hand_calculation():
    s, sites = _summary()
    h = s["hwe"]
    # Site 4: 4 hom-alt, 1 hom-ref, 0 het -> P = 5/45 = 1/9. Sites 1-3 are fully determined (p = 1).
    assert h["n_tested"] == 4 and h["n_skipped_non_diploid_or_multiallelic"] == 1
    assert h["min_p"] == pytest.approx(1 / 9)
    assert h["n_below_threshold"] == 0 and h["fraction_below_0.05"] == 0.0
    assert h["p_histogram"] == [0, 1, 0, 0, 0, 0, 0, 0, 0, 3]
    assert [round(x["hwe_p"], 4) if x["hwe_p"] is not None else None for x in sites] == [
        1.0, 1.0, 1.0, 0.1111, None]


def test_per_sample_qc_values():
    s, _ = _summary()
    q = s["sample_qc"]
    assert q["S1"]["mean_dp"] == pytest.approx(17.5) and q["S1"]["mean_gq"] == pytest.approx(52.5)
    assert q["S5"]["mean_dp"] == pytest.approx(100 / 3)  # site 1 missing -> 3 values
    assert q["S5"]["mean_gq"] == pytest.approx(200 / 3)
    assert q["S5"]["missing_rate"] == pytest.approx(0.2) and q["S1"]["missing_rate"] == 0.0
    assert [q[n]["het_rate"] for n in ("S1", "S2", "S3", "S4", "S5")] == pytest.approx(
        [0.2, 0.4, 0.2, 0.2, 0.0])
    assert [q[n]["singletons"] for n in ("S1", "S2", "S3", "S4", "S5")] == [0, 1, 0, 1, 0]
    # Ti/Tv over each sample's non-reference calls at biallelic SNPs
    assert q["S1"]["ti_tv"] == 1.0 and q["S2"]["ti_tv"] == 2.0 and q["S3"]["ti_tv"] == 1.0
    assert q["S4"]["ti_tv"] is None and q["S5"]["ti_tv"] is None


def test_outlier_flags_need_variation_and_enough_samples():
    s, _ = _summary()
    assert s["flagged_samples"] == {}  # MAD is zero here, so nothing can be called an outlier
    assert robust_z([1.0, 2.0, 3.0]) == [0.0, 0.0, 0.0]  # fewer than 5 values
    z = robust_z([10, 10.1, 9.9, 10.05, 9.95, 14.0])
    assert z[-1] > 10 and all(abs(v) < 2 for v in z[:-1])


def test_site_sink_rows():
    _, sites = _summary()
    assert len(sites) == 5
    row = sites[1]
    assert (row["chrom"], row["pos"], row["ac"], row["af_obs"]) == ("chr1", 200, 4, 0.4)
    assert sites[0]["missing_frac"] == pytest.approx(0.2)
    assert sites[4]["ac"] is None and sites[4]["hwe_p"] is None  # multiallelic


# ---------------- simulator ----------------

def test_expected_ti_tv_formula_limits():
    assert expected_ti_tv(2.1, 0.0) == pytest.approx(2.1)
    assert expected_ti_tv(2.1, 1.0) == pytest.approx(0.5)  # purely random substitutions
    values = [expected_ti_tv(2.1, f / 10) for f in range(11)]
    assert values == sorted(values, reverse=True)


def test_simulation_is_reproducible_and_truth_is_consistent():
    a, ta = simulate_vcf(n_samples=20, n_sites=300, error_fraction=0.1, excess_het_fraction=0.05,
                         n_contaminated=1, n_low_coverage=2, seed=4)
    b, tb = simulate_vcf(n_samples=20, n_sites=300, error_fraction=0.1, excess_het_fraction=0.05,
                         n_contaminated=1, n_low_coverage=2, seed=4)
    assert a == b and ta == tb
    assert (ta.n_real_sites, ta.n_random_error_sites, ta.n_excess_het_sites) == (255, 30, 15)
    assert len(ta.contaminated_samples) == 1 and len(ta.low_coverage_samples) == 2
    assert not ta.contaminated_samples & ta.low_coverage_samples
    assert sum("ARTEFACT" in l for l in a if not l.startswith("#")) == 45
    with pytest.raises(ValueError):
        simulate_vcf(n_samples=3, n_contaminated=2, n_low_coverage=2)
    with pytest.raises(ValueError):
        simulate_vcf(error_fraction=0.7, excess_het_fraction=0.4)


def test_qc_recovers_planted_truth_end_to_end():
    n_sites, error_f, excess_f = 2000, 0.2, 0.03
    lines, truth = simulate_vcf(n_samples=80, n_sites=n_sites, error_fraction=error_f,
                                excess_het_fraction=excess_f, excess_het_rate=0.9,
                                n_contaminated=2, n_low_coverage=2,
                                contamination_het_rate=0.3, seed=21)
    sites = []
    samples, variants = parse_vcf(l + "\n" for l in lines)
    s = summarize(samples, variants, site_sink=sites)

    # planted problem samples are found, and nothing else is flagged
    assert set(s["flagged_samples"]) == truth.problem_samples

    # Ti/Tv matches the analytic expectation: random-ALT sites (errors AND the
    # excess-heterozygosity artefacts) both pull it toward 0.5
    assert s["ti_tv"] == pytest.approx(expected_ti_tv(2.1, error_f + excess_f), abs=0.15)

    # HWE: the planted excess-het sites fail and essentially nothing else does.
    # (Random-error sites are carried by 1-2 samples, which is consistent with HWE.)
    failing = [x for x in sites if x["hwe_p"] is not None and x["hwe_p"] < 1e-6]
    n_excess = truth.n_excess_het_sites
    assert 0.9 * n_excess <= len(failing) <= n_excess + 2
    assert all("ARTEFACT" in x["info"] for x in failing)
