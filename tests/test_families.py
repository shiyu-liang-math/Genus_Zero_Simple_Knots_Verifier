from __future__ import annotations

from collections import Counter

import pytest

from g0sv.families import (
    family_records,
    independent_family_records,
    independent_prediction_records,
    independent_presentation_predictions,
    prediction_records,
    presentation_predictions,
)


def _identity(record):
    return (
        record.q,
        record.m,
        record.theorem_case,
        record.branch,
        record.d,
        record.epsilon,
        record.n,
        record.predicted_cs,
    )


@pytest.mark.parametrize("q", range(2, 101))
def test_divisor_generators_match_independent_scan(q: int) -> None:
    generated = Counter(map(_identity, family_records(q)))
    scanned = Counter(map(_identity, independent_family_records(q)))
    assert generated == scanned
    assert presentation_predictions(q) == independent_presentation_predictions(q)


def test_q2_boundary_includes_normalized_m1() -> None:
    records = family_records(2)
    assert any(
        record.theorem_case == 1
        and record.branch == "gcd2"
        and record.d == 2
        and record.m == 1
        for record in records
    )


def test_proper_divisor_boundaries_and_n_values() -> None:
    records = [record for q in range(2, 80) for record in family_records(q)]
    assert all(record.d != 1 for record in records if record.theorem_case in (2, 3))
    assert all(record.n != 1 for record in records if record.theorem_case in (2, 3))
    assert any(record.theorem_case == 2 and record.n == 3 for record in records)
    assert any(record.theorem_case == 3 and record.n == 2 for record in records)
    assert any(record.theorem_case == 3 and record.n and record.n > 1 and record.n % 2 for record in records)


def test_q7_overlap_is_preserved_and_unioned() -> None:
    records = family_records(7)
    by_m: dict[int, list] = {}
    for record in records:
        by_m.setdefault(record.m, []).append(record)
    overlaps = {m: group for m, group in by_m.items() if len(group) > 1}
    assert overlaps
    predictions = presentation_predictions(7)
    for m, group in overlaps.items():
        expected = set().union(*(set(record.predicted_cs) for record in group))
        assert predictions[m].issuperset(expected)

    # The case-3 exceptional pair at m=18 is transported by core exchange
    # to the missing pair +/-3 at m=30.  This is the precise place where the
    # raw README table is not yet the full marked-presentation union.
    assert predictions[18] == set(range(1, 7))
    assert predictions[30] == set(range(1, 7))
    transported = [
        record
        for record in prediction_records(7)
        if record.source_m == 18
        and record.m == 30
        and record.theorem_case == 3
        and record.n == 2
        and record.coordinate_map == "core_exchange"
    ]
    assert len(transported) == 1
    assert transported[0].c_multiplier == 2
    assert set(transported[0].predicted_cs) == set(range(1, 7))


@pytest.mark.parametrize("q", range(2, 101))
def test_prediction_transport_matches_independent_scan(q: int) -> None:
    generated = Counter(map(_identity_prediction, prediction_records(q)))
    scanned = Counter(map(_identity_prediction, independent_prediction_records(q)))
    assert generated == scanned


def _identity_prediction(record):
    return (
        record.q,
        record.m,
        record.source_m,
        record.theorem_case,
        record.branch,
        record.d,
        record.epsilon,
        record.n,
        record.coordinate_map,
        record.c_multiplier,
        record.predicted_cs,
    )
