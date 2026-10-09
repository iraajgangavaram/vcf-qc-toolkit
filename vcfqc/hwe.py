"""Exact test for Hardy-Weinberg equilibrium at a biallelic site.

Implements the recursion of Wigginton, Cutler and Abecasis (2005), which runs
in O(n) time per site and is numerically stable for large sample sizes. The
test conditions on the observed allele counts and sums the probabilities of all
heterozygote counts that are no more likely than the one observed.
"""

from __future__ import annotations

from typing import List


def hwe_exact_pvalue(n_het: int, n_hom_a: int, n_hom_b: int) -> float:
    """Two-sided exact HWE p-value from genotype counts.

    Parameters are the number of heterozygotes and of the two homozygote
    classes (in either order; the test is symmetric). Returns 1.0 for an
    empty site. Small p-values indicate too many or too few heterozygotes,
    which in a call set most often points to genotyping or mapping artefacts.
    """
    if min(n_het, n_hom_a, n_hom_b) < 0:
        raise ValueError("genotype counts must be non-negative")
    n = n_het + n_hom_a + n_hom_b
    if n == 0:
        return 1.0

    hom_rare, hom_common = sorted((n_hom_a, n_hom_b))
    rare_copies = 2 * hom_rare + n_het  # copies of the rarer allele
    if rare_copies == 0:
        return 1.0  # monomorphic: nothing to test

    probs: List[float] = [0.0] * (rare_copies + 1)

    # Start at the most probable heterozygote count (same parity as rare_copies)
    mid = int(rare_copies * (2 * n - rare_copies) / (2 * n))
    if (mid % 2) != (rare_copies % 2):
        mid += 1
    probs[mid] = 1.0
    total = 1.0

    # Walk downwards from the mode
    het, hr = mid, (rare_copies - mid) // 2
    hc = n - het - hr
    while het >= 2:
        probs[het - 2] = probs[het] * het * (het - 1.0) / (4.0 * (hr + 1.0) * (hc + 1.0))
        total += probs[het - 2]
        hr += 1
        hc += 1
        het -= 2

    # Walk upwards from the mode
    het, hr = mid, (rare_copies - mid) // 2
    hc = n - het - hr
    while het <= rare_copies - 2:
        probs[het + 2] = probs[het] * 4.0 * hr * hc / ((het + 2.0) * (het + 1.0))
        total += probs[het + 2]
        hr -= 1
        hc -= 1
        het += 2

    target = probs[n_het] * (1.0 + 1e-12)  # tolerance for floating-point ties
    p = sum(x for x in probs if x <= target)
    return min(1.0, p / total)
