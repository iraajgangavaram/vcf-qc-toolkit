"""Validate the QC metrics against simulated call sets with planted, known problems.

    python scripts/validate_qc.py

Three experiments, all seeded and reproducible:
  A. Ti/Tv versus the share of random-error sites, against the analytic expectation.
  B. HWE exact-test sensitivity to planted excess-heterozygosity sites, and its
     false-positive rate on genuine sites, as artefact strength and sample size vary.
  C. Detection of contaminated and low-coverage samples by robust outlier flags.

Writes docs/results/*.csv, docs/results/validation_summary.json and docs/figures/*.png.
"""

from __future__ import annotations

import json
import math
import pathlib
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import _style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from vcfqc import parse_vcf, summarize  # noqa: E402
from vcfqc.simulate import expected_ti_tv, simulate_vcf  # noqa: E402

FIG = ROOT / "docs" / "figures"
RES = ROOT / "docs" / "results"
TSTV = 2.1


def run(lines, **kw):
    sites: list = []
    samples, variants = parse_vcf(l + "\n" for l in lines)
    return summarize(samples, variants, site_sink=sites, **kw), sites


def exp_a_titv(seeds=range(5), n_samples=40, n_sites=6000):
    rows = []
    for f in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5):
        obs, by_bin = [], []
        for sd in seeds:
            lines, _ = simulate_vcf(n_samples=n_samples, n_sites=n_sites, tstv=TSTV,
                                    error_fraction=f, seed=100 + sd)
            s, _ = run(lines)
            obs.append(s["ti_tv"])
            by_bin.append([b["ti_tv"] for b in s["ti_tv_by_af_bin"]])
        rows.append({"error_fraction": f, "observed_mean": statistics.fmean(obs),
                     "observed_sd": statistics.stdev(obs), "expected": expected_ti_tv(TSTV, f),
                     "by_af_bin": by_bin})
    return rows


def exp_b_hwe(n_sites=3000):
    rows = []
    for n in (50, 100, 200):
        for rate in (0.5, 0.6, 0.7, 0.8, 0.9):
            lines, _ = simulate_vcf(n_samples=n, n_sites=n_sites, excess_het_fraction=0.05,
                                    excess_het_rate=rate, seed=7)
            _, sites = run(lines)
            art = [x for x in sites if "ARTEFACT" in x["info"] and x["hwe_p"] is not None]
            real = [x for x in sites if "ARTEFACT" not in x["info"] and x["hwe_p"] is not None]
            rows.append({
                "n_samples": n, "artefact_het_rate": rate,
                "n_artefact_sites": len(art), "n_real_sites": len(real),
                "sensitivity_p<1e-6": sum(x["hwe_p"] < 1e-6 for x in art) / len(art),
                "sensitivity_p<1e-3": sum(x["hwe_p"] < 1e-3 for x in art) / len(art),
                "false_positive_p<1e-6": sum(x["hwe_p"] < 1e-6 for x in real) / len(real),
                "false_positive_p<0.05": sum(x["hwe_p"] < 0.05 for x in real) / len(real),
            })
    return rows


def qq_data(excess_fraction, n_samples=100, n_sites=4000, seed=3):
    lines, _ = simulate_vcf(n_samples=n_samples, n_sites=n_sites,
                            excess_het_fraction=excess_fraction, seed=seed)
    _, sites = run(lines)
    ps = sorted(x["hwe_p"] for x in sites if x["hwe_p"] is not None)
    n = len(ps)
    exp = [-math.log10((i + 0.5) / n) for i in range(n)]
    obs = [-math.log10(max(p, 1e-300)) for p in ps]
    return exp[::-1], obs[::-1]


def exp_c_samples(seeds=range(10), n_samples=100, n_sites=2000):
    rows = []
    for rate in (0.02, 0.05, 0.1, 0.2, 0.3):
        det_c = tot_c = det_l = tot_l = false_flags = clean = 0
        for sd in seeds:
            lines, truth = simulate_vcf(n_samples=n_samples, n_sites=n_sites, n_contaminated=2,
                                        n_low_coverage=2, contamination_het_rate=rate, seed=500 + sd)
            s, _ = run(lines)
            flagged = set(s["flagged_samples"])
            det_c += len(flagged & truth.contaminated_samples)
            tot_c += len(truth.contaminated_samples)
            det_l += len(flagged & truth.low_coverage_samples)
            tot_l += len(truth.low_coverage_samples)
            false_flags += len(flagged - truth.problem_samples)
            clean += n_samples - len(truth.problem_samples)
        rows.append({"contamination_het_rate": rate,
                     "contaminated_detected": det_c / tot_c,
                     "low_coverage_detected": det_l / tot_l,
                     "clean_samples_falsely_flagged": false_flags / clean,
                     "n_runs": len(list(seeds))})
    return rows


def scatter_data():
    lines, truth = simulate_vcf(n_samples=100, n_sites=2000, n_contaminated=2, n_low_coverage=2,
                                contamination_het_rate=0.3, seed=11)
    s, _ = run(lines)
    pts = []
    for name, q in s["sample_qc"].items():
        kind = ("contaminated" if name in truth.contaminated_samples
                else "low coverage" if name in truth.low_coverage_samples else "clean")
        pts.append((q["missing_rate"], q["het_rate"], kind, name in s["flagged_samples"]))
    return pts


def main():
    _style.apply()
    FIG.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)

    a = exp_a_titv()
    b = exp_b_hwe()
    c = exp_c_samples()

    import csv
    def write(name, rows, drop=()):
        keys = [k for k in rows[0] if k not in drop]
        with open(RES / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items() if k in keys})
    write("validation_titv.csv", a, drop=("by_af_bin",))
    write("validation_hwe.csv", b)
    write("validation_samples.csv", c)

    # ---- Figure 1: Ti/Tv ----
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    xs = [r["error_fraction"] * 100 for r in a]
    ax.plot(xs, [r["expected"] for r in a], color=_style.INK_2, linestyle="--", marker="s",
            label="Expected (analytic)", markeredgecolor=_style.SURFACE, markeredgewidth=1.2)
    ax.errorbar(xs, [r["observed_mean"] for r in a], yerr=[r["observed_sd"] for r in a],
                color=_style.BLUE, marker="o", capsize=3, label="Observed (mean +/- SD, 5 runs)",
                markeredgecolor=_style.SURFACE, markeredgewidth=1.2)
    ax.axhline(0.5, color=_style.GRID, linewidth=1)
    ax.annotate("random errors: 0.5", (xs[-1], 0.5), xytext=(-4, 6), textcoords="offset points",
                ha="right", color=_style.INK_2, fontsize=8.5)
    ax.set_title("Ti/Tv falls as false positives enter")
    ax.set_xlabel("Random-error sites (% of call set)")
    ax.set_ylabel("Ti/Tv")
    ax.set_ylim(0, 2.5)
    ax.legend(loc="upper right")

    ax = axes[1]
    centers = [(i + 0.5) / 10 for i in range(10)]
    for r, color, marker, label in ((a[0], _style.BLUE, "o", "0% errors"),
                                    (a[3], _style.ORANGE, "s", "20% errors")):
        means = []
        for i in range(10):
            vals = [row[i] for row in r["by_af_bin"] if row[i] is not None]
            means.append(statistics.fmean(vals) if vals else float("nan"))
        ax.plot(centers, means, color=color, marker=marker, label=label,
                markeredgecolor=_style.SURFACE, markeredgewidth=1.2)
    ax.set_title("Errors concentrate at low allele frequency")
    ax.set_xlabel("Alternate allele frequency (observed)")
    ax.set_ylabel("Ti/Tv")
    ax.set_ylim(0, 3.2)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "titv_validation.png", dpi=200)
    plt.close(fig)

    # ---- Figure 2: HWE ----
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    for frac, color, marker, label in ((0.0, _style.BLUE, "o", "Clean call set"),
                                       (0.05, _style.ORANGE, "s", "5% excess-het artefacts")):
        e, o = qq_data(frac)
        ax.scatter(e, o, s=9, color=color, marker=marker, label=label, linewidths=0, alpha=0.8)
    lim = ax.get_xlim()[1]
    ax.plot([0, lim], [0, lim], color=_style.INK_2, linestyle="--", linewidth=1)
    ax.set_title("HWE p-values: clean data follow the diagonal")
    ax.set_xlabel("Expected -log10 p")
    ax.set_ylabel("Observed -log10 p")
    ax.legend(loc="upper left")

    ax = axes[1]
    for n, color, marker in ((50, _style.BLUE, "o"), (100, _style.ORANGE, "s"), (200, _style.AQUA, "^")):
        rows = [r for r in b if r["n_samples"] == n]
        ax.plot([r["artefact_het_rate"] for r in rows], [r["sensitivity_p<1e-6"] for r in rows],
                color=color, marker=marker, label=f"{n} samples",
                markeredgecolor=_style.SURFACE, markeredgewidth=1.2)
    ax.set_title("Detection depends on artefact strength")
    ax.set_xlabel("Heterozygous share at artefact sites")
    ax.set_ylabel("Share of artefact sites with p < 1e-6")
    ax.set_ylim(-0.03, 1.08)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "hwe_validation.png", dpi=200)
    plt.close(fig)

    # ---- Figure 3: sample QC ----
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    style = {"clean": (_style.BLUE, "o", 22), "contaminated": (_style.ORANGE, "s", 60),
             "low coverage": (_style.AQUA, "^", 60)}
    for kind, (color, marker, size) in style.items():
        pts = [p for p in scatter_data() if p[2] == kind]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=size, color=color, marker=marker,
                   label=f"{kind} ({len(pts)})", edgecolors=_style.SURFACE, linewidths=0.8,
                   alpha=0.85 if kind == "clean" else 1.0)
    ax.set_title("Problem samples separate on two axes")
    ax.set_xlabel("Missing-genotype rate")
    ax.set_ylabel("Heterozygosity rate")
    ax.legend(loc="center right")

    ax = axes[1]
    xs = [r["contamination_het_rate"] for r in c]
    ax.plot(xs, [r["contaminated_detected"] for r in c], color=_style.ORANGE, marker="s",
            label="Contaminated samples detected", markeredgecolor=_style.SURFACE, markeredgewidth=1.2)
    ax.plot(xs, [r["low_coverage_detected"] for r in c], color=_style.AQUA, marker="^",
            label="Low-coverage samples detected", markeredgecolor=_style.SURFACE, markeredgewidth=1.2)
    ax.plot(xs, [r["clean_samples_falsely_flagged"] for r in c], color=_style.BLUE, marker="o",
            label="Clean samples falsely flagged", markeredgecolor=_style.SURFACE, markeredgewidth=1.2)
    ax.set_title("Detection versus contamination level")
    ax.set_xlabel("Extra heterozygous calls in contaminated samples")
    ax.set_ylabel("Share of samples")
    ax.set_ylim(-0.04, 1.12)
    ax.legend(loc="center right")
    fig.tight_layout()
    fig.savefig(FIG / "sample_qc_validation.png", dpi=200)
    plt.close(fig)

    summary = {
        "titv_max_abs_deviation_from_expected": max(abs(r["observed_mean"] - r["expected"]) for r in a),
        "hwe_false_positive_rate_p<1e-6_max": max(r["false_positive_p<1e-6"] for r in b),
        "hwe_false_positive_rate_p<0.05_max": max(r["false_positive_p<0.05"] for r in b),
    }
    (RES / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    for name, rows in (("A", a), ("B", b), ("C", c)):
        print("---", name)
        for r in rows:
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k != "by_af_bin"})


if __name__ == "__main__":
    main()
