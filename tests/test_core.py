from __future__ import annotations

import json
from math import gcd
from pathlib import Path

import pytest

from g0sv.core import optimized_extrema, reference_walk, reference_walk_points, surface_genus


def test_readme_l25_6_example() -> None:
    expected = reference_walk(5, 6, 1)
    assert expected.min_g == -4
    assert expected.max_g == 4
    assert expected.argmin == 6
    assert expected.argmax == 4
    assert expected.width == 8
    assert optimized_extrema(5, 6, 1) == expected
    assert surface_genus(expected.width, 5) == 0


def test_stored_l25_6_example() -> None:
    path = Path(__file__).resolve().parents[1] / "examples" / "l25_6.json"
    stored = json.loads(path.read_text())
    extrema = reference_walk(stored["q"], stored["m"], stored["c"])
    assert (extrema.min_g, extrema.max_g, extrema.argmin, extrema.argmax, extrema.width) == (
        stored["min_g"],
        stored["max_g"],
        stored["argmin"],
        stored["argmax"],
        stored["width"],
    )


@pytest.mark.parametrize("q", range(2, 16))
def test_all_small_units_against_reference(q: int) -> None:
    p = q * q
    for m in range(1, p):
        if gcd(m, p) != 1:
            continue
        for c in range(1, q):
            if gcd(c, q) == 1:
                assert optimized_extrema(q, m, c) == reference_walk(q, m, c), (q, m, c)


def test_walk_wraparound_and_endpoint() -> None:
    points = list(reference_walk_points(7, 48, 6))
    assert points[0].x == 0 and points[0].g == 0
    assert points[-1].x == 49 and points[-1].g == 0
    assert len(points) == 50
    extrema = reference_walk(7, 48, 6)
    assert extrema.argmin < 49
    assert extrema.argmax < 49


def test_both_signs_have_expected_width() -> None:
    for c in (1, 4):
        extrema = optimized_extrema(5, 6, c)
        assert extrema.width == 8
        assert surface_genus(extrema.width, 5) == 0


def test_tied_extrema_choose_earliest_position() -> None:
    for q, m, c in ((3, 1, 1), (4, 5, 3), (8, 17, 1), (9, 37, 8)):
        assert optimized_extrema(q, m, c) == reference_walk(q, m, c)


def test_surface_genus_rejects_invalid_width() -> None:
    with pytest.raises(ValueError):
        surface_genus(7, 5)
    with pytest.raises(ValueError):
        surface_genus(4, 5)
