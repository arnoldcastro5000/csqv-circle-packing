"""Tests for tools/validate_pck.py, the independent AS-WRITTEN .pck submission gate.

These lock that the gate passes a real sealed record, catches a geometric overlap or wall
violation, and rejects a file that is non-overlapping but Packomania-noncompliant (unsorted
radii, a mis-declared max, a non-positive radius). main() must return a nonzero exit code on
any invalid file so `validate && submit` is safe.
"""

from __future__ import annotations

import importlib.util
import os
import types

import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_SEALED_N98 = os.path.join(_ROOT, "discoveries", "2026-09-15-22-csqv-98", "csqv98.pck")
_SEALED_N98_SUM = 5.208554826179


def _load_tool(name: str) -> types.ModuleType:
    path = os.path.join(_ROOT, "tools", f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_pck(path: str, rows: list[tuple[float, float, float]]) -> None:
    largest = max(r for _, _, r in rows)
    with open(path, "w") as f:
        f.write(f"{largest:.15f}\ntest\n")
        for x, y, r in rows:
            f.write(f"{x:.15f} {y:.15f} {r:.15f}\n")


def test_passes_sealed_record() -> None:
    vp = _load_tool("validate_pck")
    res = vp.validate(_SEALED_N98)
    assert res["feasible"] is True
    assert res["valid"] is True
    assert res["sum_radii"] == pytest.approx(_SEALED_N98_SUM, abs=1e-11)


def test_catches_overlap(tmp_path) -> None:
    vp = _load_tool("validate_pck")
    path = os.path.join(tmp_path, "csqv2.pck")
    _write_pck(path, [(-0.05, 0.0, 0.06), (0.05, 0.0, 0.06)])  # 0.02 overlap
    res = vp.validate(path)
    assert res["feasible"] is False
    assert res["valid"] is False
    assert res["worst_overlap"] == pytest.approx(0.02, abs=1e-12)


def test_catches_wall_violation(tmp_path) -> None:
    vp = _load_tool("validate_pck")
    path = os.path.join(tmp_path, "csqv1.pck")
    _write_pck(path, [(0.40, 0.0, 0.15)])  # |x|+r = 0.55 > 0.5
    res = vp.validate(path)
    assert res["feasible"] is False
    assert res["worst_wall"] == pytest.approx(0.05, abs=1e-12)


def test_unsorted_is_invalid(tmp_path) -> None:
    vp = _load_tool("validate_pck")
    path = os.path.join(tmp_path, "csqv2.pck")
    _write_pck(path, [(-0.3, 0.0, 0.05), (0.3, 0.0, 0.04)])  # feasible but descending
    res = vp.validate(path)
    assert res["feasible"] is True
    assert res["sorted_incr"] is False
    assert res["valid"] is False


def test_negative_radius_is_invalid(tmp_path) -> None:
    vp = _load_tool("validate_pck")
    path = os.path.join(tmp_path, "csqv1.pck")
    _write_pck(path, [(0.0, 0.0, -0.01)])
    res = vp.validate(path)
    assert res["all_r_pos"] is False
    assert res["valid"] is False


def test_declared_max_mismatch_is_invalid(tmp_path) -> None:
    vp = _load_tool("validate_pck")
    path = os.path.join(tmp_path, "csqv1.pck")
    with open(path, "w") as f:
        f.write("0.999999999999999\ntest\n0.0 0.0 0.05\n")
    res = vp.validate(path)
    assert res["declared_max_diff"] > 1e-6
    assert res["valid"] is False


def test_main_exit_code(tmp_path) -> None:
    vp = _load_tool("validate_pck")
    good = os.path.join(tmp_path, "csqv2.pck")
    _write_pck(good, [(-0.3, 0.0, 0.04), (0.3, 0.0, 0.05)])
    bad = os.path.join(tmp_path, "csqv2bad.pck")
    _write_pck(bad, [(-0.05, 0.0, 0.06), (0.05, 0.0, 0.06)])
    assert vp.main([good]) == 0
    assert vp.main([bad]) == 1


def test_short_header_raises(tmp_path) -> None:
    vp = _load_tool("validate_pck")
    path = os.path.join(tmp_path, "csqv0.pck")
    with open(path, "w") as f:
        f.write("0.05\n")
    with pytest.raises(ValueError):
        vp.load_pck(path)
