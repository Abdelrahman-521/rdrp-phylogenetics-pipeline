"""Unit tests for the pipeline's Python steps. No network needed: NCBI calls
are replaced with fakes."""

import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"
DATA = Path(__file__).resolve().parents[1] / "data"
sys.path.insert(0, str(BIN))

import compare_trees  # noqa: E402
import extract_domain  # noqa: E402
import phylobot_utils as u  # noqa: E402
import rename_headers  # noqa: E402
from repair_fasta import repair  # noqa: E402


# ---------- parsing ----------

def test_parse_fasta_drops_ncbi_error_text():
    text = (
        ">CAA25029.1 polyprotein [Cowpea mosaic virus]\nMKLV\nAGLK\n"
        "Error: F a i l e d  t o  u n d e r s t a n d  i d:  MG833031.1\n"
        ">BAA03151.1 polyprotein\nVVGR\n"
    )
    records, junk = u.parse_fasta_text(text)
    assert [u.record_id(h) for h, _ in records] == ["CAA25029.1", "BAA03151.1"]
    assert records[0][1] == "MKLVAGLK"
    assert len(junk) == 1 and junk[0].startswith("Error")


def test_empty_header_is_junk_not_a_record():
    records, junk = u.parse_fasta_text(">lonely header\n>real\nMK\n")
    assert [h for h, _ in records] == ["real"]
    assert junk == [">lonely header"]


@pytest.mark.parametrize("rec_id, expected", [
    ("CAA25029.1", "CAA25029"),
    ("BAA21898.1", "BAA21898"),
    ("lcl|MG833031.1_prot_AWC26954.1_1", "MG833031"),
])
def test_source_accession(rec_id, expected):
    assert u.source_accession(rec_id) == expected


def test_missing_accessions_handles_versions_and_cds_records():
    records = [("CAA25029.1 x", "MK"), ("BAA21898.1 y", "MK"), ("lcl|LC310707.1_prot_BBC14926.1_1 z", "MK")]
    assert u.missing_accessions(records, ["CAA25029.1", "BAA21898", "LC310707", "KU645789"]) == ["KU645789"]


# ---------- repair (the checkerrors.py logic) ----------

def test_repair_recovers_from_nuccore_and_reports():
    raw = ">CAA25029.1 a\nMKLV\n+ some curl error line\n"
    calls = []

    def fake_fetch(acc):
        calls.append(acc)
        return f">lcl|{acc}.1_prot_XYZ.1_1 [gene=pol]\nMSTR\n" if acc == "MG833031.1" else ""

    records, report = repair(raw, ["CAA25029.1", "MG833031.1", "NOPE123"], fake_fetch)
    assert calls == ["MG833031.1", "NOPE123"]
    assert report["recovered"] == ["MG833031.1"]
    assert report["still_missing"] == ["NOPE123"]
    assert report["junk_lines"] == 1
    assert len(records) == 2


def test_repair_does_nothing_when_complete():
    records, report = repair(">A1.1\nMK\n", ["A1.1"], lambda acc: pytest.fail("should not fetch"))
    assert report == {"junk_lines": 0, "recovered": [], "still_missing": []}


# ---------- domain extraction (the get_superfamily.py logic) ----------

def test_extract_cuts_one_based_inclusive_coordinates(tmp_path):
    regions_file = tmp_path / "r.tsv"
    regions_file.write_text("seq_id\tsignature\tstart\tend\nP1.1\tSSF56672\t3\t6\nP2.1\tOTHER\t1\t2\n")
    regions = extract_domain.load_regions(regions_file, "regions", "SSF56672")
    out, problems = extract_domain.extract([("P1.1 desc", "ABCDEFGH"), ("P2.1", "XYZ")], regions, "SSF56672")
    assert out == [("P1.1_3_6 SSF56672", "CDEF")]
    assert problems == []


def test_extract_reads_raw_interproscan_tsv(tmp_path):
    tsv = tmp_path / "ips.tsv"
    tsv.write_text("P1.1\tmd5\t100\tSUPERFAMILY\tSSF56672\tDNA/RNA polymerases\t2\t4\t1e-50\n")
    regions = extract_domain.load_regions(tsv, "interproscan", "SSF56672")
    assert dict(regions) == {"P1.1": [(2, 4)]}


def test_extract_flags_regions_past_the_end_and_missing_records():
    regions = {"P1.1": [(5, 50)], "GONE.1": [(1, 2)]}
    out, problems = extract_domain.extract([("P1.1", "ABCDEFGH")], regions, "SSF56672")
    assert out == []
    assert any("past the end" in p for p in problems)
    assert any("GONE.1" in p for p in problems)


def test_bundled_regions_cover_every_accession():
    regions = extract_domain.load_regions(DATA / "rdrp_domain_regions.tsv", "regions", "SSF56672")
    accessions = u.read_list(DATA / "accessions.txt")
    covered = {u.source_accession(rid) for rid in regions}
    assert {u.accession_base(a) for a in accessions} == covered


# ---------- renaming ----------

def test_rename_uses_short_names_and_removes_spaces():
    names = {"CAA25029.1": "CPMV", "YP_009336629.1": "HuAV 1"}
    out, unnamed = rename_headers.rename(
        [("CAA25029.1_1178_1670 SSF56672", "MK"), ("YP_009336629.1_2300_2816 SSF56672", "MK"), ("X9.1_1_2", "MK")],
        names,
    )
    assert [h for h, _ in out] == ["CPMV", "HuAV_1", "X9.1_1_2"]
    assert unnamed == ["X9.1_1_2"]


def test_rename_refuses_duplicate_names():
    with pytest.raises(ValueError):
        rename_headers.rename([("A.1_1_2", "M"), ("B.1_1_2", "M")], {"A.1": "Same", "B.1": "Same"})


def test_bundled_names_cover_bundled_domains():
    names = {r["seq_id"]: r["short_name"] for r in u.read_tsv(DATA / "virus_names.tsv")}
    records, _ = u.parse_fasta_text((DATA / "SSF56672_sequences.fasta").read_text())
    _, unnamed = rename_headers.rename(records, names)
    assert unnamed == []


# ---------- tree comparison ----------

def _tree(newick):
    from io import StringIO
    from Bio import Phylo
    return Phylo.read(StringIO(newick), "newick")


def test_identical_trees_have_rf_zero_even_when_rooted_differently():
    a = _tree("((A,B)90,(C,D)80,(E,F)70);")
    b = _tree("(A,(B,((C,D)80,(E,F)70)));")
    result = compare_trees.compare(a, b)
    assert result["rf"] == 0
    assert result["shared_splits"] == 3


def test_rf_counts_conflicting_splits_and_support():
    a = _tree("((A,B)100,(C,D)60,(E,F)90);")
    b = _tree("((A,C)55,(B,D)50,(E,F)95);")
    result = compare_trees.compare(a, b)
    assert result["rf"] == 4          # AB and CD vs AC and BD
    assert result["rf_normalised"] == pytest.approx(4 / 6)
    assert result["support_shared_a"] == 90
    assert result["support_only_a"] == pytest.approx(80)


def test_changed_relatives():
    a = _tree("(((A:1,B:1):1,C:2):1,(D:1,E:1):2);")
    b = _tree("(((A:1,C:1):1,B:2):1,(D:1,E:1):2);")
    moved = compare_trees.changed_relatives(a, b)
    assert set(moved) == {"A", "B", "C"}


def test_alignment_stats(tmp_path):
    aln = tmp_path / "a.fasta"
    aln.write_text(">x\nAC-D\n>y\nACED\n")
    stats = compare_trees.alignment_stats(aln)
    assert stats["columns"] == 4
    assert stats["gap_fraction"] == pytest.approx(1 / 8)
    assert stats["mean_pairwise_identity"] == pytest.approx(1.0)
