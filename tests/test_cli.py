import gzip
import json

from vcfqc.cli import main

VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    "chr1\t10\t.\tA\tG\t50\tPASS\tAF=0.4\tGT\t0/1\n"
    "chr1\t20\t.\tC\tA\t30\tPASS\tAF=0.6\tGT\t1/1\n"
)


def test_text_report_plain_and_gzipped(tmp_path, capsys):
    plain = tmp_path / "x.vcf"
    plain.write_text(VCF)
    gz = tmp_path / "x.vcf.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(VCF)

    reports = []
    for path in (plain, gz):
        assert main([str(path)]) == 0
        reports.append(capsys.readouterr().out)
    assert reports[0] == reports[1]
    assert "Records:            2" in reports[0]
    assert "Ti/Tv (SNPs):       1.000" in reports[0]


def test_json_output(tmp_path, capsys):
    path = tmp_path / "x.vcf"
    path.write_text(VCF)
    assert main([str(path), "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["n_records"] == 2 and data["samples"]["S1"]["het"] == 1


def test_errors_return_nonzero_not_traceback(tmp_path, capsys):
    assert main([str(tmp_path / "missing.vcf")]) == 1
    bad = tmp_path / "bad.vcf"
    bad.write_text("not a vcf\n")
    assert main([str(bad)]) == 1
    assert "error" in capsys.readouterr().err


def test_extended_report_and_site_table(tmp_path, capsys):
    from vcfqc.simulate import simulate_vcf

    lines, truth = simulate_vcf(n_samples=40, n_sites=400, n_contaminated=1, n_low_coverage=1,
                                contamination_het_rate=0.4, seed=2)
    vcf = tmp_path / "sim.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    table = tmp_path / "sites.tsv"

    assert main([str(vcf), "--site-table", str(table), "--hwe-threshold", "1e-4"]) == 0
    out = capsys.readouterr().out
    assert "HWE (exact test):" in out and "p < 0.0001" in out
    assert "Flagged samples (robust outliers):" in out
    for name in truth.problem_samples:
        assert name in out.split("Flagged samples (robust outliers):")[1]

    rows = table.read_text().strip().split("\n")
    assert rows[0] == "chrom\tpos\tref\talt\tqual\tmissing_frac\tac\taf_obs\thwe_p"
    assert len(rows) == 401  # header + one row per site


def test_json_contains_new_sections(tmp_path, capsys):
    path = tmp_path / "x.vcf"
    path.write_text(VCF)
    assert main([str(path), "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    for key in ("hwe", "site_qc", "sample_qc", "substitution_spectrum", "ti_tv_by_af_bin"):
        assert key in data
    assert data["hwe"]["n_tested"] == 2
