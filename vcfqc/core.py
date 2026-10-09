"""Parse VCF text and compute standard variant-call QC summaries."""

from __future__ import annotations

import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CT")
BASES = frozenset("ACGT")


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alts: List[str]
    qual: Optional[float]
    filter: str
    info: Dict[str, str]
    genotypes: List[str] = field(default_factory=list)  # raw GT string per sample


def parse_vcf(lines: Iterable[str]) -> Tuple[List[str], Iterator[Variant]]:
    """Return (sample_names, variant_iterator). Consumes header lines eagerly.

    Raises ValueError if the '#CHROM' header line is missing or a data line is
    malformed (fewer than 8 columns, non-integer POS).
    """
    it = iter(lines)
    samples: List[str] = []
    for line in it:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            cols = line.rstrip("\n").split("\t")
            samples = cols[9:]
            break
        raise ValueError("VCF is missing the '#CHROM' header line")
    else:
        raise ValueError("VCF is empty or has no '#CHROM' header line")

    def variants() -> Iterator[Variant]:
        for lineno, line in enumerate(it, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) < 8:
                raise ValueError(f"data line {lineno} has {len(f)} columns, expected >= 8")
            try:
                pos = int(f[1])
            except ValueError:
                raise ValueError(f"data line {lineno}: POS {f[1]!r} is not an integer") from None
            info = {}
            if f[7] != ".":
                for item in f[7].split(";"):
                    k, _, v = item.partition("=")
                    info[k] = v
            gts = []
            if len(f) > 9 and "GT" in f[8].split(":"):
                gi = f[8].split(":").index("GT")
                gts = [s.split(":")[gi] if len(s.split(":")) > gi else "." for s in f[9:]]
            yield Variant(
                chrom=f[0],
                pos=pos,
                ref=f[3].upper(),
                alts=f[4].upper().split(",") if f[4] != "." else [],
                qual=None if f[5] == "." else float(f[5]),
                filter=f[6],
                info=info,
                genotypes=gts,
            )

    return samples, variants()


def classify(ref: str, alt: str) -> str:
    """Classify one REF/ALT pair: SNP, MNP, INS, DEL, SYMBOLIC or OTHER."""
    if alt.startswith("<") or alt == "*" or "[" in alt or "]" in alt:
        return "SYMBOLIC"
    if not (set(ref) <= BASES and set(alt) <= BASES):
        return "OTHER"
    if len(ref) == len(alt):
        return "SNP" if len(ref) == 1 else "MNP"
    return "INS" if len(alt) > len(ref) else "DEL"


def is_transition(ref: str, alt: str) -> bool:
    """True for A<->G or C<->T. Only meaningful for single-base substitutions."""
    return (ref in PURINES and alt in PURINES) or (ref in PYRIMIDINES and alt in PYRIMIDINES)


def genotype_class(gt: str) -> str:
    """'missing', 'hom_ref', 'hom_alt' or 'het' for a GT string such as 0/1 or 1|2."""
    alleles = gt.replace("|", "/").split("/")
    if any(a in (".", "") for a in alleles):
        return "missing"
    if len(set(alleles)) > 1:
        return "het"
    return "hom_ref" if alleles[0] == "0" else "hom_alt"


def summarize(samples: List[str], variants: Iterable[Variant], af_bins: int = 10) -> dict:
    """Aggregate QC statistics over all variant records."""
    n_records = 0
    types: Counter = Counter()
    per_chrom: Counter = Counter()
    filters: Counter = Counter()
    quals: List[float] = []
    indel_lengths: Counter = Counter()
    ti = tv = 0
    af_hist = [0] * af_bins
    sample_gt: Dict[str, Counter] = defaultdict(Counter)

    for v in variants:
        n_records += 1
        per_chrom[v.chrom] += 1
        filters["PASS" if v.filter in ("PASS", ".") else "FAIL"] += 1
        if v.qual is not None:
            quals.append(v.qual)
        for alt in v.alts:
            kind = classify(v.ref, alt)
            types[kind] += 1
            if kind == "SNP":
                if is_transition(v.ref, alt):
                    ti += 1
                else:
                    tv += 1
            elif kind in ("INS", "DEL"):
                indel_lengths[len(alt) - len(v.ref)] += 1
        raw_af = v.info.get("AF", "").split(",")[0]
        try:
            af = float(raw_af)
            if 0.0 <= af <= 1.0:
                af_hist[min(int(af * af_bins + 1e-9), af_bins - 1)] += 1
        except ValueError:
            pass
        for name, gt in zip(samples, v.genotypes):
            sample_gt[name][genotype_class(gt)] += 1

    per_sample = {}
    for name in samples:
        c = sample_gt[name]
        per_sample[name] = {
            "hom_ref": c["hom_ref"], "het": c["het"],
            "hom_alt": c["hom_alt"], "missing": c["missing"],
            "het_hom_alt_ratio": c["het"] / c["hom_alt"] if c["hom_alt"] else None,
        }

    return {
        "n_records": n_records,
        "variant_types": dict(types),
        "ti_tv": ti / tv if tv else None,
        "n_transitions": ti,
        "n_transversions": tv,
        "per_chromosome": dict(per_chrom),
        "filter": dict(filters),
        "qual_mean": statistics.fmean(quals) if quals else None,
        "qual_median": statistics.median(quals) if quals else None,
        "indel_length_distribution": dict(sorted(indel_lengths.items())),
        "af_histogram": af_hist,
        "samples": per_sample,
    }
