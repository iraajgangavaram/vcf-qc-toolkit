"""Command-line interface: python -m vcfqc FILE.vcf[.gz] [--json] [--site-table sites.tsv]."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys

from .core import parse_vcf, summarize


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _fmt(x, nd=3):
    return "NA" if x is None else f"{x:.{nd}f}"


def format_report(s: dict) -> str:
    out = [
        f"Records:            {s['n_records']}",
        "Variant types:      " + ", ".join(f"{k}={v}" for k, v in sorted(s["variant_types"].items())),
        f"Ti/Tv (SNPs):       {_fmt(s['ti_tv'])}  ({s['n_transitions']} ts / {s['n_transversions']} tv)",
        "FILTER:             " + ", ".join(f"{k}={v}" for k, v in sorted(s["filter"].items())),
        f"QUAL mean/median:   {_fmt(s['qual_mean'], 2)} / {_fmt(s['qual_median'], 2)}",
        "Per chromosome:     " + ", ".join(f"{k}={v}" for k, v in s["per_chromosome"].items()),
        "Indel lengths:      "
        + (", ".join(f"{k:+d}:{v}" for k, v in s["indel_length_distribution"].items()) or "none"),
        "AF histogram:       " + " ".join(str(x) for x in s["af_histogram"]),
        "Substitutions:      " + ", ".join(f"{k}={v}" for k, v in s["substitution_spectrum"].items()),
    ]
    if s["samples"]:
        site, hwe = s["site_qc"], s["hwe"]
        out += [
            "",
            f"Site QC:            mean missingness {_fmt(site['mean_site_missingness'], 4)}, "
            f"{site['n_sites_missingness_above_threshold']} sites above "
            f"{site['missingness_threshold']:.0%}, singletons {site['singleton_sites']}, "
            f"doubletons {site['doubleton_sites']}",
            f"HWE (exact test):   {hwe['n_tested']} sites tested, "
            f"{hwe['n_below_threshold']} with p < {hwe['threshold']:g}, "
            f"min p {('NA' if hwe['min_p'] is None else format(hwe['min_p'], '.2e'))}, "
            f"{hwe['n_skipped_non_diploid_or_multiallelic']} skipped",
            "Observed AF hist:   " + " ".join(str(x) for x in s["af_histogram_observed"]),
            "Ti/Tv by AF bin:    " + " ".join(_fmt(b["ti_tv"], 2) for b in s["ti_tv_by_af_bin"]),
            "",
            "Sample              hom_ref   het  hom_alt  missing  het/hom_alt",
        ]
        for name, c in s["samples"].items():
            out.append(
                f"{name:<18}{c['hom_ref']:>8}{c['het']:>6}{c['hom_alt']:>9}{c['missing']:>9}"
                f"{_fmt(c['het_hom_alt_ratio'], 2):>13}"
            )
        out += ["", "Sample QC           missing   het_rate  mean_DP  mean_GQ  singletons  Ti/Tv"]
        for name, q in s["sample_qc"].items():
            out.append(
                f"{name:<18}{_fmt(q['missing_rate'], 4):>9}{_fmt(q['het_rate'], 4):>11}"
                f"{_fmt(q['mean_dp'], 1):>9}{_fmt(q['mean_gq'], 1):>9}{q['singletons']:>12}"
                f"{_fmt(q['ti_tv'], 2):>7}"
            )
        if s["flagged_samples"]:
            out += ["", "Flagged samples (robust outliers):"]
            out += [f"  {n}: {'; '.join(r)}" for n, r in s["flagged_samples"].items()]
        else:
            out += ["", "Flagged samples:    none"]
    return "\n".join(out)


def write_site_table(path: str, sites: list) -> None:
    cols = ["chrom", "pos", "ref", "alt", "qual", "missing_frac", "ac", "af_obs", "hwe_p"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for row in sites:
            w.writerow(["NA" if row[c] is None else row[c] for c in cols])


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="vcfqc", description=__doc__)
    p.add_argument("vcf")
    p.add_argument("--json", action="store_true", help="emit JSON instead of a text report")
    p.add_argument("--site-table", metavar="TSV",
                   help="also write one row per site (missingness, AC, AF, HWE p-value)")
    p.add_argument("--hwe-threshold", type=float, default=1e-6,
                   help="HWE p-value below which a site counts as failing (default 1e-6)")
    p.add_argument("--missingness-threshold", type=float, default=0.1,
                   help="site missingness above which a site is counted (default 0.1)")
    args = p.parse_args(argv)
    sites: list = []
    try:
        with _open(args.vcf) as fh:
            samples, variants = parse_vcf(fh)
            summary = summarize(
                samples, variants, hwe_threshold=args.hwe_threshold,
                missingness_threshold=args.missingness_threshold,
                site_sink=sites if args.site_table else None,
            )
        if args.site_table:
            write_site_table(args.site_table, sites)
    except (OSError, ValueError) as e:
        print(f"vcfqc: error: {e}", file=sys.stderr)
        return 1
    print(json.dumps(summary, indent=2) if args.json else format_report(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
