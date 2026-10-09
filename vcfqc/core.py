"""Parse VCF text and compute variant-call quality-control summaries.

Two levels of output are produced in a single streaming pass:

* call-set level: variant types, Ti/Tv, FILTER/QUAL, allele-frequency spectra,
  substitution spectrum, Ti/Tv by allele frequency, Hardy-Weinberg tests;
* sample level: genotype counts, missingness, heterozygosity, depth, genotype
  quality, singletons, and robust outlier flags.
"""

from __future__ import annotations

import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

from .hwe import hwe_exact_pvalue

PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CT")
BASES = frozenset("ACGT")
_COMP = str.maketrans("ACGT", "TGCA")


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
    fmt: List[str] = field(default_factory=list)  # FORMAT keys, e.g. ['GT', 'DP', 'GQ']
    sample_data: List[List[str]] = field(default_factory=list)  # per-sample FORMAT values

    def sample_value(self, sample_index: int, key: str) -> Optional[float]:
        """Numeric FORMAT value (e.g. DP, GQ) for one sample, or None if absent/missing."""
        if key not in self.fmt or sample_index >= len(self.sample_data):
            return None
        fields = self.sample_data[sample_index]
        i = self.fmt.index(key)
        if i >= len(fields) or fields[i] in (".", ""):
            return None
        try:
            return float(fields[i])
        except ValueError:
            return None


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
            fmt: List[str] = []
            sample_data: List[List[str]] = []
            gts: List[str] = []
            if len(f) > 9:
                fmt = f[8].split(":")
                sample_data = [s.split(":") for s in f[9:]]
                if "GT" in fmt:
                    gi = fmt.index("GT")
                    gts = [d[gi] if len(d) > gi else "." for d in sample_data]
            yield Variant(
                chrom=f[0],
                pos=pos,
                ref=f[3].upper(),
                alts=f[4].upper().split(",") if f[4] != "." else [],
                qual=None if f[5] == "." else float(f[5]),
                filter=f[6],
                info=info,
                genotypes=gts,
                fmt=fmt,
                sample_data=sample_data,
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


def substitution_class(ref: str, alt: str) -> str:
    """Strand-collapsed substitution class: one of C>A, C>G, C>T, T>A, T>C, T>G.

    A purine reference (A or G) is reverse-complemented so that, for example,
    G>A and C>T are counted together.
    """
    if ref in PURINES:
        ref, alt = ref.translate(_COMP), alt.translate(_COMP)
    return f"{ref}>{alt}"


SUBSTITUTION_CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")


def genotype_class(gt: str) -> str:
    """'missing', 'hom_ref', 'hom_alt' or 'het' for a GT string such as 0/1 or 1|2."""
    alleles = gt.replace("|", "/").split("/")
    if any(a in (".", "") for a in alleles):
        return "missing"
    if len(set(alleles)) > 1:
        return "het"
    return "hom_ref" if alleles[0] == "0" else "hom_alt"


def _is_diploid(gt: str) -> bool:
    return len(gt.replace("|", "/").split("/")) == 2


def robust_z(values: List[float]) -> List[float]:
    """Robust z-scores using the median and MAD (scaled to match a normal SD).

    Returns zeros if the MAD is zero or fewer than 5 values are given, because
    outliers cannot be assessed reliably in either case.
    """
    if len(values) < 5:
        return [0.0] * len(values)
    med = statistics.median(values)
    mad = statistics.median(abs(v - med) for v in values) * 1.4826
    if mad == 0:
        return [0.0] * len(values)
    return [(v - med) / mad for v in values]


def summarize(
    samples: List[str],
    variants: Iterable[Variant],
    af_bins: int = 10,
    hwe_threshold: float = 1e-6,
    missingness_threshold: float = 0.1,
    outlier_z: float = 3.5,
    site_sink: Optional[List[dict]] = None,
) -> dict:
    """Aggregate QC statistics over all variant records in one pass.

    Parameters
    ----------
    hwe_threshold : sites with an exact HWE p-value below this are counted as failing.
    missingness_threshold : sites with a larger fraction of missing genotypes are counted.
    outlier_z : samples with a robust |z| above this on heterozygosity, missingness or
        depth are flagged.
    site_sink : if a list is given, one dict per variant record with genotypes is
        appended (chrom, pos, ref, alt, qual, missing_frac, ac, af_obs, hwe_p, info),
        so sites can be inspected or filtered individually.
    """
    n_records = 0
    types: Counter = Counter()
    per_chrom: Counter = Counter()
    filters: Counter = Counter()
    quals: List[float] = []
    indel_lengths: Counter = Counter()
    ti = tv = 0
    af_hist = [0] * af_bins
    af_hist_obs = [0] * af_bins
    sample_gt: Dict[str, Counter] = defaultdict(Counter)
    spectrum: Counter = Counter()
    bin_ts = [0] * af_bins
    bin_tv = [0] * af_bins

    n_sites_gt = 0
    site_missing_sum = 0.0
    n_high_missing = 0
    singleton_sites = doubleton_sites = 0
    hwe_p: List[float] = []
    hwe_skipped = 0
    s_ts: Dict[str, int] = defaultdict(int)
    s_tv: Dict[str, int] = defaultdict(int)
    s_singletons: Dict[str, int] = defaultdict(int)
    dp_sum: Dict[str, float] = defaultdict(float)
    dp_n: Dict[str, int] = defaultdict(int)
    gq_sum: Dict[str, float] = defaultdict(float)
    gq_n: Dict[str, int] = defaultdict(int)

    def bin_of(af: float) -> int:
        return min(int(af * af_bins + 1e-9), af_bins - 1)

    for v in variants:
        n_records += 1
        per_chrom[v.chrom] += 1
        filters["PASS" if v.filter in ("PASS", ".") else "FAIL"] += 1
        if v.qual is not None:
            quals.append(v.qual)

        info_af: Optional[float] = None
        try:
            cand = float(v.info.get("AF", "").split(",")[0])
            if 0.0 <= cand <= 1.0:
                info_af = cand
                af_hist[bin_of(cand)] += 1
        except ValueError:
            pass

        # ---- genotype-level, site-wise statistics ----
        gt_af: Optional[float] = None
        site_ac: Optional[int] = None
        site_hwe_p: Optional[float] = None
        miss_frac = None
        classes: List[str] = []
        biallelic = len(v.alts) == 1
        if samples and v.genotypes:
            n_sites_gt += 1
            classes = [genotype_class(gt) for gt in v.genotypes[: len(samples)]]
            for name, cls in zip(samples, classes):
                sample_gt[name][cls] += 1
            counts = Counter(classes)
            miss_frac = counts["missing"] / len(samples)
            site_missing_sum += miss_frac
            if miss_frac > missingness_threshold:
                n_high_missing += 1

            diploid = all(_is_diploid(g) for g, c in zip(v.genotypes, classes) if c != "missing")
            alleles_ok = all(
                set(g.replace("|", "/").split("/")) <= {"0", "1"}
                for g, c in zip(v.genotypes, classes) if c != "missing"
            )
            called = len(samples) - counts["missing"]
            if biallelic and diploid and alleles_ok and called > 0:
                ac = 2 * counts["hom_alt"] + counts["het"]
                gt_af = ac / (2 * called)
                af_hist_obs[bin_of(gt_af)] += 1
                site_ac = ac
                site_hwe_p = hwe_exact_pvalue(counts["het"], counts["hom_ref"], counts["hom_alt"])
                hwe_p.append(site_hwe_p)
                if ac == 1:
                    singleton_sites += 1
                    s_singletons[samples[classes.index("het")]] += 1
                elif ac == 2:
                    doubleton_sites += 1
            else:
                hwe_skipped += 1  # multiallelic, non-diploid, or no called genotypes

            for i, name in enumerate(samples):
                dp = v.sample_value(i, "DP")
                if dp is not None:
                    dp_sum[name] += dp
                    dp_n[name] += 1
                gq = v.sample_value(i, "GQ")
                if gq is not None:
                    gq_sum[name] += gq
                    gq_n[name] += 1

        if site_sink is not None and classes:
            site_sink.append({
                "chrom": v.chrom, "pos": v.pos, "ref": v.ref, "alt": ",".join(v.alts),
                "qual": v.qual, "missing_frac": miss_frac, "ac": site_ac,
                "af_obs": gt_af, "hwe_p": site_hwe_p, "info": v.info,
            })

        # ---- per-alternate-allele statistics ----
        for alt in v.alts:
            kind = classify(v.ref, alt)
            types[kind] += 1
            if kind == "SNP":
                transition = is_transition(v.ref, alt)
                ti += transition
                tv += not transition
                spectrum[substitution_class(v.ref, alt)] += 1
                af_for_bin = gt_af if gt_af is not None else info_af
                if biallelic and af_for_bin is not None:
                    if transition:
                        bin_ts[bin_of(af_for_bin)] += 1
                    else:
                        bin_tv[bin_of(af_for_bin)] += 1
                if biallelic and classes:
                    for name, cls in zip(samples, classes):
                        if cls in ("het", "hom_alt"):
                            if transition:
                                s_ts[name] += 1
                            else:
                                s_tv[name] += 1
            elif kind in ("INS", "DEL"):
                indel_lengths[len(alt) - len(v.ref)] += 1

    per_sample = {}
    sample_qc = {}
    for name in samples:
        c = sample_gt[name]
        per_sample[name] = {
            "hom_ref": c["hom_ref"], "het": c["het"],
            "hom_alt": c["hom_alt"], "missing": c["missing"],
            "het_hom_alt_ratio": c["het"] / c["hom_alt"] if c["hom_alt"] else None,
        }
        total = sum(c.values())
        called = total - c["missing"]
        sample_qc[name] = {
            "missing_rate": c["missing"] / total if total else None,
            "het_rate": c["het"] / called if called else None,
            "mean_dp": dp_sum[name] / dp_n[name] if dp_n[name] else None,
            "mean_gq": gq_sum[name] / gq_n[name] if gq_n[name] else None,
            "singletons": s_singletons[name],
            "ti_tv": s_ts[name] / s_tv[name] if s_tv[name] else None,
        }

    flagged: Dict[str, List[str]] = {}
    for metric in ("het_rate", "missing_rate", "mean_dp"):
        usable = [n for n in samples if sample_qc[n][metric] is not None]
        zs = robust_z([sample_qc[n][metric] for n in usable])
        for n, z in zip(usable, zs):
            sample_qc[n][f"{metric}_z"] = z
            if abs(z) > outlier_z:
                flagged.setdefault(n, []).append(f"{metric} z={z:+.1f}")

    hwe_hist = [0] * 10
    for p in hwe_p:
        hwe_hist[min(int(p * 10), 9)] += 1
    hwe = {
        "n_tested": len(hwe_p),
        "n_skipped_non_diploid_or_multiallelic": hwe_skipped,
        "threshold": hwe_threshold,
        "n_below_threshold": sum(p < hwe_threshold for p in hwe_p),
        "fraction_below_0.05": (sum(p < 0.05 for p in hwe_p) / len(hwe_p)) if hwe_p else None,
        "min_p": min(hwe_p) if hwe_p else None,
        "p_histogram": hwe_hist,
    }

    ti_tv_by_bin = []
    for b in range(af_bins):
        ti_tv_by_bin.append({
            "af_low": b / af_bins, "af_high": (b + 1) / af_bins,
            "transitions": bin_ts[b], "transversions": bin_tv[b],
            "ti_tv": bin_ts[b] / bin_tv[b] if bin_tv[b] else None,
        })

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
        "substitution_spectrum": {k: spectrum[k] for k in SUBSTITUTION_CLASSES},
        "af_histogram_observed": af_hist_obs,
        "ti_tv_by_af_bin": ti_tv_by_bin,
        "site_qc": {
            "n_sites_with_genotypes": n_sites_gt,
            "mean_site_missingness": site_missing_sum / n_sites_gt if n_sites_gt else None,
            "n_sites_missingness_above_threshold": n_high_missing,
            "missingness_threshold": missingness_threshold,
            "singleton_sites": singleton_sites,
            "doubleton_sites": doubleton_sites,
        },
        "hwe": hwe,
        "sample_qc": sample_qc,
        "flagged_samples": flagged,
    }
