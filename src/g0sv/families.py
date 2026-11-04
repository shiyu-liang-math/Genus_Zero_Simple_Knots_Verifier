"""The four positive-presentation generators and BBL residue predictions."""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import Iterable, Iterator


@dataclass(frozen=True, slots=True)
class FamilyRecord:
    q: int
    m: int
    theorem_case: int
    branch: str
    d: int
    epsilon: int | None
    n: int | None
    predicted_cs: tuple[int, ...]

    @property
    def family(self) -> str:
        return f"case{self.theorem_case}_{self.branch}"


@dataclass(frozen=True, slots=True)
class PredictionRecord:
    """A family contribution transported into a marked presentation."""

    q: int
    m: int
    source_m: int
    theorem_case: int
    branch: str
    d: int
    epsilon: int | None
    n: int | None
    coordinate_map: str
    c_multiplier: int
    predicted_cs: tuple[int, ...]

    @property
    def family(self) -> str:
        return f"case{self.theorem_case}_{self.branch}"


def signed_residues(q: int, values: Iterable[int]) -> tuple[int, ...]:
    residues: set[int] = set()
    for value in values:
        positive = value % q
        negative = (-value) % q
        if positive:
            residues.add(positive)
        if negative:
            residues.add(negative)
    result = tuple(sorted(residues))
    if any(gcd(c, q) != 1 for c in result):
        raise AssertionError(f"BBL prediction contains a non-order-q residue for q={q}: {result}")
    return result


def case1_gcd1(q: int) -> Iterator[FamilyRecord]:
    p = q * q
    for d in range(1, q + 1):
        if gcd(d, q) == 1:
            yield FamilyRecord(q, (d * q + 1) % p, 1, "gcd1", d, None, None, signed_residues(q, (1, d)))


def case1_gcd2(q: int) -> Iterator[FamilyRecord]:
    p = q * q
    for d in range(1, q + 1):
        if gcd(d, q) == 2:
            yield FamilyRecord(q, (d * q + 1) % p, 1, "gcd2", d, None, None, signed_residues(q, (1,)))


def _odd_divisors(number: int) -> Iterator[int]:
    for divisor in range(3, number + 1, 2):
        if number % divisor == 0:
            yield divisor


def case2(q: int) -> Iterator[FamilyRecord]:
    p = q * q
    for epsilon in (-1, 1):
        number = 2 * q + epsilon
        for d in _odd_divisors(number):
            n = number // d
            if n >= 3 and n % 2 == 1:
                m = d * (q - epsilon) % p
                yield FamilyRecord(q, m, 2, "odd_2q", d, epsilon, n, signed_residues(q, (1, d)))


def case3(q: int) -> Iterator[FamilyRecord]:
    p = q * q
    for epsilon in (-1, 1):
        number = q - epsilon
        for d in _odd_divisors(number):
            n = number // d
            if n > 1:
                values = (1, d, 2) if n == 2 else (1, d)
                m = d * (q - epsilon) % p
                yield FamilyRecord(q, m, 3, "odd_q", d, epsilon, n, signed_residues(q, values))


def family_records(q: int) -> list[FamilyRecord]:
    if q < 2:
        return []
    return [*case1_gcd1(q), *case1_gcd2(q), *case2(q), *case3(q)]


def prediction_records(q: int) -> list[PredictionRecord]:
    """Transport every paper class into every generated equivalent marking.

    Besides identical marked presentations, the paper uses the mirror
    ``m -> -m`` and core exchange ``m -> m^{-1}``.  A mirror fixes ``mu`` up
    to the already-retained sign, while a core exchange sends source ``mu``
    to target ``mu' = target_m * mu``.  The q=7 overlap between m=18 and
    m=30 is the smallest place where this closure adds classes.
    """

    raw = family_records(q)
    targets = sorted({record.m for record in raw})
    p = q * q
    result: list[PredictionRecord] = []
    for record in raw:
        inverse = pow(record.m, -1, p)
        for target in targets:
            if target == record.m:
                coordinate_map = "identity"
                multiplier = 1
            elif target == (-record.m) % p:
                coordinate_map = "mirror"
                multiplier = 1
            elif target == inverse:
                coordinate_map = "core_exchange"
                multiplier = target % q
            elif target == (-inverse) % p:
                coordinate_map = "core_exchange_mirror"
                multiplier = target % q
            else:
                continue
            transported = signed_residues(q, (multiplier * c for c in record.predicted_cs))
            result.append(
                PredictionRecord(
                    q=q,
                    m=target,
                    source_m=record.m,
                    theorem_case=record.theorem_case,
                    branch=record.branch,
                    d=record.d,
                    epsilon=record.epsilon,
                    n=record.n,
                    coordinate_map=coordinate_map,
                    c_multiplier=multiplier,
                    predicted_cs=transported,
                )
            )
    return result


def presentation_predictions(q: int) -> dict[int, set[int]]:
    predictions: dict[int, set[int]] = {}
    for record in prediction_records(q):
        predictions.setdefault(record.m, set()).update(record.predicted_cs)
    return predictions


def units(q: int) -> list[int]:
    return [c for c in range(1, q) if gcd(c, q) == 1]


def independent_family_records(q: int) -> list[FamilyRecord]:
    """Scan all candidate d and test the defining equations directly.

    This deliberately does not call the divisor-based production generators.
    It is used for coverage verification and generator boundary tests.
    """

    p = q * q
    records: list[FamilyRecord] = []
    for d in range(1, q + 1):
        common = gcd(d, q)
        if common == 1:
            records.append(FamilyRecord(q, (d * q + 1) % p, 1, "gcd1", d, None, None, signed_residues(q, (1, d))))
        if common == 2:
            records.append(FamilyRecord(q, (d * q + 1) % p, 1, "gcd2", d, None, None, signed_residues(q, (1,))))

    for epsilon in (-1, 1):
        for d in range(3, 2 * q + abs(epsilon) + 1):
            if d % 2 == 0:
                continue
            number = 2 * q + epsilon
            if number == d * (number // d):
                n = number // d
                if n >= 3 and n % 2 == 1:
                    records.append(FamilyRecord(q, d * (q - epsilon) % p, 2, "odd_2q", d, epsilon, n, signed_residues(q, (1, d))))
        for d in range(3, q + abs(epsilon) + 1):
            if d % 2 == 0:
                continue
            number = q - epsilon
            if number == d * (number // d):
                n = number // d
                if n > 1:
                    values = (1, d, 2) if n == 2 else (1, d)
                    records.append(FamilyRecord(q, d * (q - epsilon) % p, 3, "odd_q", d, epsilon, n, signed_residues(q, values)))
    return records


def independent_prediction_records(q: int) -> list[PredictionRecord]:
    """Independently scan equivalent targets and induced residue maps."""

    raw = independent_family_records(q)
    targets = sorted({record.m for record in raw})
    p = q * q
    result: list[PredictionRecord] = []
    for record in raw:
        inverse = pow(record.m, -1, p)
        for target in targets:
            if target == record.m:
                coordinate_map, multiplier = "identity", 1
            elif (target + record.m) % p == 0:
                coordinate_map, multiplier = "mirror", 1
            elif (target * record.m) % p == 1:
                coordinate_map, multiplier = "core_exchange", target % q
            elif (target * record.m) % p == p - 1:
                coordinate_map, multiplier = "core_exchange_mirror", target % q
            else:
                continue
            transported = signed_residues(q, [multiplier * c for c in record.predicted_cs])
            result.append(
                PredictionRecord(
                    q,
                    target,
                    record.m,
                    record.theorem_case,
                    record.branch,
                    record.d,
                    record.epsilon,
                    record.n,
                    coordinate_map,
                    multiplier,
                    transported,
                )
            )
    return result


def independent_presentation_predictions(q: int) -> dict[int, set[int]]:
    predictions: dict[int, set[int]] = {}
    for record in independent_prediction_records(q):
        predictions.setdefault(record.m, set()).update(record.predicted_cs)
    return predictions
