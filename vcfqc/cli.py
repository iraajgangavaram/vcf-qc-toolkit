"""Command-line interface: python -m vcfqc FILE.vcf[.gz] [--json]."""

from __future__ import annotations

import argparse
import gzip
import json
import sys

from .core import parse_vcf, summarize


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def format_report(s: dict) -> str:
    fmt = lambda x, nd=3: "NA" if x is None else f"{x:.{nd}f}"
    out = [
        f"Records:            {s['n_records']}",
        f"Variant types:      " + ", ".join(f"{k}={v}" for k, v in sorted(s["variant_types"].items())),
        f"Ti/Tv (SNPs):       {fmt(s['ti_tv'])}  ({s['n_transitions']} ts / {s['n_transversions']} tv)",
        f"FILTER:             " + ", ".join(f"{k}={v}" for k, v in sorted(s["filter"].items())),
        f"QUAL mean/median:   {fmt(s['qual_mean'], 2)} / {fmt(s['qual_median'], 2)}",
        "Per chromosome:     " + ", ".join(f"{k}={v}" for k, v in s["per_chromosome"].items()),
        "Indel lengths:      " + (", ".join(f"{k:+d}:{v}" for k, v in s["indel_length_distribution"].items()) or "none"),
        "AF histogram:       " + " ".join(str(x) for x in s["af_histogram"]),
    ]
    if s["samples"]:
        out.append("\nSample              hom_ref   het  hom_alt  missing  het/hom_alt")
        for name, c in s["samples"].items():
            out.append(
                f"{name:<18}{c['hom_ref']:>8}{c['het']:>6}{c['hom_alt']:>9}{c['missing']:>9}"
                f"{fmt(c['het_hom_alt_ratio'], 2):>13}"
            )
    return "\n".join(out)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="vcfqc", description=__doc__)
    p.add_argument("vcf")
    p.add_argument("--json", action="store_true", help="emit JSON instead of a text report")
    args = p.parse_args(argv)
    try:
        with _open(args.vcf) as fh:
            samples, variants = parse_vcf(fh)
            summary = summarize(samples, variants)
    except (OSError, ValueError) as e:
        print(f"vcfqc: error: {e}", file=sys.stderr)
        return 1
    print(json.dumps(summary, indent=2) if args.json else format_report(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
