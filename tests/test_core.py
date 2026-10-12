"""Tests for the ChimeraX-independent helpers.  Run with:  python tests/test_core.py"""
import importlib.util
import math
import os

_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "core.py")
_spec = importlib.util.spec_from_file_location("compare_core", _path)
core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(core)


def test_rmsd():
    assert core.rmsd([3.0, 4.0]) == math.sqrt(12.5)
    assert math.isnan(core.rmsd([]))


def test_tm_score_identical_is_one():
    assert abs(core.tm_score([0.0] * 100, 100) - 1.0) < 1e-12


def test_tm_score_drops_with_distance_and_coverage():
    assert core.tm_score([5.0] * 100, 100) < core.tm_score([1.0] * 100, 100)
    assert abs(core.tm_score([0.0] * 50, 100) - 0.5) < 1e-12


def test_identity_and_ss():
    assert core.percent_identity("ACDE", "ACDF") == 75.0
    assert core.ss_agreement("HHEE", "HH--") == 50.0
    assert core.ss_composition("HHE-") == (50.0, 25.0, 25.0)


def test_gapped_columns():
    cols = core.gapped_columns(4, 3, [(0, 0), (2, 1), (3, 2)])
    assert cols == [(0, 0), (1, None), (2, 1), (3, 2)]
    cols = core.gapped_columns(2, 4, [(0, 1), (1, 2)])
    assert cols == [(None, 0), (0, 1), (1, 2), (None, 3)]


def test_format_alignment():
    text = core.format_alignment("wt", "mut", "ACDE", "HHHH", "AGDE", "HH-H",
        [(0, 0), (1, 1), (2, 2), (3, 3)])
    lines = text.splitlines()
    assert lines[1].endswith("ACDE") and lines[3].endswith("AGDE")
    assert lines[2].endswith("| ||")
    assert lines[5].endswith("  ^ ")


def test_substitutions():
    subs = core.substitutions("ACDE", [10, 11, 12, 13], "AGDE", [1, 2, 3, 4],
        [(0, 0), (1, 1), (2, 2), (3, 3)], [0.1, 0.9, 0.2, 0.3])
    assert subs == [("C", 11, "G", 2, 0.9)]


def test_star_columns_and_rows():
    # ref ACDE; b lacks ref position 1 and has one extra residue before ref 3;
    # c matches ref fully.
    cols = core.star_columns(4, [(4, [(0, 0), (2, 1), (3, 3)]),
        (4, [(0, 0), (1, 1), (2, 2), (3, 3)])])
    assert cols == [(0, [0, 0]), (1, [None, 1]), (2, [1, 2]), (None, [2, None]), (3, [3, 3])]
    seq_rows, ss_rows, same_res, same_ss = core.multiple_rows(
        ["ACDE", "ADXE", "ACDF"], ["HHHH", "HH-H", "HHHE"], cols)
    assert seq_rows == ["ACD-E", "A-DXE", "ACD-F"]
    assert same_res == "* *  "
    assert same_ss == "= =  "
    assert core.multiple_stats(cols, same_res, same_ss) == (3, 2, 2)
    text = core.format_multiple(["ref", "b", "c"], seq_rows, ss_rows, same_res, same_ss)
    assert "ACD-E" in text and "identical" in text


def test_star_columns_trailing_insert():
    cols = core.star_columns(2, [(4, [(0, 1), (1, 2)])])
    assert cols == [(None, [0]), (0, [1]), (1, [2]), (None, [3])]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
