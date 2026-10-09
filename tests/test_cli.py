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
