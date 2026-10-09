"""Simulate a multi-sample VCF with known ground truth, for validating QC metrics.

Real variants
    Allele frequencies follow a 1/x site-frequency spectrum (many rare variants,
    few common ones), genotypes are drawn under Hardy-Weinberg equilibrium, and
    substitutions are transitions with probability tstv / (1 + tstv).
Artefacts (planted, known)
    * random-error sites: ALT is a uniformly random different base, so the
      expected Ti/Tv is 0.5; carried by one or two samples; low QUAL.
    * excess-heterozygosity sites: a share ``excess_het_rate`` (default 85%) of
      samples is heterozygous regardless of allele frequency, as happens when
      reads from paralogous loci collapse onto one reference position.
Problem samples (planted, known)
    * contaminated: extra heterozygous calls at real sites;
    * low coverage: lower depth and many missing genotypes.

The simulated INFO field contains an ``ARTEFACT`` flag marking planted artefacts.
It exists only so that experiments can score the QC output; ``vcfqc`` never uses it.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import List, Set, Tuple

_PARTNER = {"A": "G", "G": "A", "C": "T", "T": "C"}
_TRANSVERSIONS = {"A": "CT", "G": "CT", "C": "AG", "T": "AG"}


@dataclass
class Truth:
    n_real_sites: int
    n_random_error_sites: int
    n_excess_het_sites: int
    contaminated_samples: Set[str] = field(default_factory=set)
    low_coverage_samples: Set[str] = field(default_factory=set)

    @property
    def problem_samples(self) -> Set[str]:
        return self.contaminated_samples | self.low_coverage_samples


def expected_ti_tv(tstv: float, error_fraction: float) -> float:
    """Analytic Ti/Tv of a call set in which `error_fraction` of SNPs are random errors."""
    p_ts = tstv / (1.0 + tstv)
    f = error_fraction
    ts = (1 - f) * p_ts + f / 3.0
    tv = (1 - f) * (1 - p_ts) + f * 2.0 / 3.0
    return ts / tv


def _sample_af(rng: random.Random, lo: float = 0.005, hi: float = 0.95) -> float:
    """Draw from a density proportional to 1/x on [lo, hi]."""
    return lo * (hi / lo) ** rng.random()


def simulate_vcf(
    n_samples: int = 60,
    n_sites: int = 5000,
    tstv: float = 2.1,
    error_fraction: float = 0.0,
    excess_het_fraction: float = 0.0,
    n_contaminated: int = 0,
    n_low_coverage: int = 0,
    contamination_het_rate: float = 0.3,
    excess_het_rate: float = 0.85,
    mean_depth: float = 30.0,
    seed: int = 0,
) -> Tuple[List[str], Truth]:
    """Return (VCF lines, ground truth).

    ``error_fraction`` and ``excess_het_fraction`` are shares of all sites.
    """
    if n_contaminated + n_low_coverage > n_samples:
        raise ValueError("more problem samples than samples")
    if error_fraction + excess_het_fraction >= 1.0:
        raise ValueError("artefact fractions must sum to less than 1")
    rng = random.Random(seed)
    names = [f"S{i + 1:03d}" for i in range(n_samples)]
    order = list(range(n_samples))
    rng.shuffle(order)
    contaminated = {names[i] for i in order[:n_contaminated]}
    low_cov = {names[i] for i in order[n_contaminated:n_contaminated + n_low_coverage]}

    n_err = round(error_fraction * n_sites)
    n_exc = round(excess_het_fraction * n_sites)
    kinds = ["real"] * (n_sites - n_err - n_exc) + ["error"] * n_err + ["excess"] * n_exc
    rng.shuffle(kinds)

    p_ts = tstv / (1.0 + tstv)
    chroms = ["chr1", "chr2", "chr3"]
    lines = [
        "##fileformat=VCFv4.2",
        "##source=vcfqc.simulate",
        *[f"##contig=<ID={c}>" for c in chroms],
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##INFO=<ID=ARTEFACT,Number=0,Type=Flag,Description="Planted artefact (simulation truth)">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *names]),
    ]

    depth_scale = {n: (0.3 if n in low_cov else rng.uniform(0.8, 1.2)) for n in names}
    miss_rate = {n: (0.25 if n in low_cov else 0.01) for n in names}
    pos = 0
    for site_index, kind in enumerate(kinds):
        pos += rng.randint(50, 500)
        chrom = chroms[min(site_index * len(chroms) // n_sites, len(chroms) - 1)]
        ref = rng.choice("ACGT")
        n_alt_alleles = 0
        gts: List[int] = []  # 0, 1, 2 copies of ALT; -1 = decide later

        if kind == "real":
            if rng.random() < p_ts:
                alt = _PARTNER[ref]
            else:
                alt = rng.choice(_TRANSVERSIONS[ref])
            q = _sample_af(rng)
            for n in names:
                if n in contaminated and rng.random() < contamination_het_rate:
                    gts.append(1)
                else:
                    gts.append((rng.random() < q) + (rng.random() < q))
            qual = max(10.0, rng.gauss(70, 20))
        elif kind == "error":
            alt = rng.choice([b for b in "ACGT" if b != ref])
            gts = [0] * n_samples
            for i in rng.sample(range(n_samples), 1 if rng.random() < 0.7 else 2):
                gts[i] = 1
            qual = max(5.0, rng.gauss(25, 10))
        else:  # excess heterozygosity
            alt = rng.choice([b for b in "ACGT" if b != ref])
            gts = [1 if rng.random() < excess_het_rate else rng.choice((0, 2)) for _ in names]
            qual = max(10.0, rng.gauss(50, 15))

        fields = []
        alt_count = called = 0
        for gi, n in enumerate(names):
            if rng.random() < miss_rate[n]:
                fields.append("./.:.:.")
                continue
            g = gts[gi]
            alt_count += g
            called += 2
            gt = {0: "0/0", 1: "0/1" if rng.random() < 0.5 else "1/0", 2: "1/1"}[g]
            dp = max(1, round(rng.gauss(mean_depth * depth_scale[n], math.sqrt(mean_depth))))
            gq = min(99, max(1, round(rng.gauss(60 if dp > 8 else 25, 15))))
            fields.append(f"{gt}:{dp}:{gq}")
        af = alt_count / called if called else 0.0
        info = f"AF={af:.4f}" + (";ARTEFACT" if kind != "real" else "")
        lines.append("\t".join([chrom, str(pos), ".", ref, alt, f"{qual:.1f}", "PASS", info,
                                "GT:DP:GQ", *fields]))

    return lines, Truth(
        n_real_sites=n_sites - n_err - n_exc,
        n_random_error_sites=n_err,
        n_excess_het_sites=n_exc,
        contaminated_samples=contaminated,
        low_coverage_samples=low_cov,
    )
