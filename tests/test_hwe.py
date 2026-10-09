from fractions import Fraction
from math import factorial

import pytest

from vcfqc.hwe import hwe_exact_pvalue


def _brute_force_p(n_het, n_hom_a, n_hom_b):
    """Exact HWE p-value by enumerating every genotype table with the same allele
    counts and using exact rational arithmetic (no recursion, no floating point).

    P(n_AB | n, n_A) = n! * n_A! * n_B! * 2^n_AB / (n_AA! n_AB! n_BB! (2n)!)
    """
    n = n_het + n_hom_a + n_hom_b
    n_a = 2 * n_hom_a + n_het
    n_b = 2 * n_hom_b + n_het

    def prob(h):
        if (n_a - h) % 2 or (n_b - h) % 2 or h > min(n_a, n_b):
            return None
        aa, bb = (n_a - h) // 2, (n_b - h) // 2
        if aa < 0 or bb < 0:
            return None
        return Fraction(
            factorial(n) * factorial(n_a) * factorial(n_b) * 2**h,
            factorial(aa) * factorial(h) * factorial(bb) * factorial(2 * n),
        )

    table = {h: prob(h) for h in range(0, min(n_a, n_b) + 1)}
    table = {h: p for h, p in table.items() if p is not None}
    assert sum(table.values()) == 1  # sanity check on the enumeration itself
    observed = table[n_het]
    return float(sum(p for p in table.values() if p <= observed))


def test_matches_exact_enumeration_for_every_table_up_to_n_25():
    checked = 0
    for n in range(1, 26):
        for hom_a in range(n + 1):
            for het in range(n - hom_a + 1):
                hom_b = n - hom_a - het
                got = hwe_exact_pvalue(het, hom_a, hom_b)
                assert got == pytest.approx(_brute_force_p(het, hom_a, hom_b), abs=1e-10), (
                    het, hom_a, hom_b)
                checked += 1
    assert checked > 3000


def test_symmetric_in_the_two_homozygote_classes():
    for het, a, b in [(10, 40, 5), (3, 90, 7), (50, 1, 49)]:
        assert hwe_exact_pvalue(het, a, b) == pytest.approx(hwe_exact_pvalue(het, b, a), abs=1e-14)


def test_perfect_equilibrium_gives_p_of_one():
    # p = q = 0.5, n = 100: expected 25 / 50 / 25
    assert hwe_exact_pvalue(50, 25, 25) == pytest.approx(1.0)


def test_excess_and_deficit_of_heterozygotes_are_both_flagged():
    assert hwe_exact_pvalue(100, 0, 0) < 1e-20  # everyone heterozygous
    assert hwe_exact_pvalue(0, 50, 50) < 1e-20  # no heterozygotes at all


def test_monomorphic_and_empty_sites_are_not_flagged():
    assert hwe_exact_pvalue(0, 0, 80) == 1.0
    assert hwe_exact_pvalue(0, 0, 0) == 1.0


def test_stable_for_biobank_scale_sample_sizes():
    # n = 300k, allele frequency 0.5 -> expected genotype counts 75k / 150k / 75k
    assert hwe_exact_pvalue(150_000, 75_000, 75_000) == pytest.approx(1.0)
    excess = hwe_exact_pvalue(250_000, 25_000, 25_000)  # far too many heterozygotes
    deficit = hwe_exact_pvalue(50_000, 125_000, 125_000)  # far too few
    assert 0.0 <= excess < 1e-10 and 0.0 <= deficit < 1e-10


def test_pvalues_are_uniform_for_data_simulated_under_hwe():
    import random

    rng = random.Random(1)
    ps = []
    for _ in range(4000):
        n, q = 200, rng.uniform(0.1, 0.5)
        het = hom_b = hom_a = 0
        for _i in range(n):
            g = (rng.random() < q) + (rng.random() < q)
            if g == 1:
                het += 1
            elif g == 2:
                hom_b += 1
            else:
                hom_a += 1
        ps.append(hwe_exact_pvalue(het, hom_a, hom_b))
    # The exact test is conservative (discrete), so rejections must not exceed nominal.
    assert sum(p < 0.05 for p in ps) / len(ps) <= 0.05
    assert sum(p < 0.01 for p in ps) / len(ps) <= 0.01


def test_rejects_negative_counts():
    with pytest.raises(ValueError):
        hwe_exact_pvalue(-1, 5, 5)
