"""Exact jump-function implementations.

The reference implementation follows Section 10.1 of the paper literally.
The production implementation instead uses the inverse-membership criterion
and block summaries, so its arithmetic and control flow are independent of
the forward-multiple construction.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import Iterator, NamedTuple


@dataclass(frozen=True, slots=True)
class Extrema:
    min_g: int
    max_g: int
    argmin: int
    argmax: int

    @property
    def width(self) -> int:
        return self.max_g - self.min_g


class WalkPoint(NamedTuple):
    x: int
    g: int
    increment: int
    in_s: bool


def least_positive_residue(value: int, modulus: int) -> int:
    """Return the representative in ``{1, ..., modulus}``."""

    residue = value % modulus
    return residue if residue else modulus


def knot_parameters(q: int, m: int, c: int) -> tuple[int, int, int, int]:
    """Return ``(p, m_star, ell, t)`` for the order-q class ``cq``."""

    if q < 2:
        raise ValueError("q must be at least 2")
    if not 1 <= c < q or gcd(c, q) != 1:
        raise ValueError("c must be a unit with 1 <= c < q")
    p = q * q
    m %= p
    if not m or gcd(m, p) != 1:
        raise ValueError("m must be a unit modulo q^2")
    m_star = pow(m, -1, p)
    ell = least_positive_residue(m_star * c * q, p)
    if ell % q:
        raise AssertionError("ell is not divisible by q for an order-q class")
    t = ell // q
    if not 0 < t < q or gcd(t, q) != 1:
        raise AssertionError("invalid order-q jump parameter")
    return p, m_star, ell, t


def reference_walk_points(q: int, m: int, c: int) -> Iterator[WalkPoint]:
    """Stream the literal normalized walk, including x=0 and x=p.

    The jump set is built from forward multiples of ``m`` exactly as in the
    README and paper, not by the inverse-membership shortcut.
    """

    p, _m_star, ell, t = knot_parameters(q, m, c)
    m %= p
    in_s = bytearray(p + 1)
    position = 0
    for _ in range(ell):
        position = (position + m) % p
        if position == 0:
            position = p
        if in_s[position]:
            raise AssertionError("forward multiples repeated before ell")
        in_s[position] = 1
        if position == p:
            position = 0
    if sum(in_s) != ell:
        raise AssertionError("jump-set cardinality mismatch")

    g_value = 0
    yield WalkPoint(0, 0, 0, False)
    for x in range(1, p + 1):
        member = bool(in_s[x])
        increment = t - q if member else t
        g_value += increment
        yield WalkPoint(x, g_value, increment, member)
    if g_value != 0:
        raise AssertionError("g(p) != g(0)")


def reference_walk(q: int, m: int, c: int) -> Extrema:
    """Compute exact extrema with the literal forward-multiple walk."""

    p = q * q
    min_g = max_g = 0
    argmin = argmax = 0
    last: WalkPoint | None = None
    for point in reference_walk_points(q, m, c):
        last = point
        if point.x == p:
            continue
        if point.g < min_g:
            min_g, argmin = point.g, point.x
        if point.g > max_g:
            max_g, argmax = point.g, point.x
    if last is None or last.x != p or last.g != 0:
        raise AssertionError("incomplete reference walk")
    return Extrema(min_g, max_g, argmin, argmax)


class _SummaryTree:
    """Segment tree of sums and earliest prefix extrema."""

    def __init__(self, values: list[int]):
        self.n = len(values)
        size = 1
        while size < self.n:
            size *= 2
        self.size = size
        count = 2 * size
        self.total = [0] * count
        self.min_prefix = [0] * count
        self.max_prefix = [0] * count
        self.argmin = [0] * count
        self.argmax = [0] * count
        self.length = [0] * count
        for index in range(size):
            node = size + index
            if index < self.n:
                self._set_leaf(node, values[index])
            else:
                self.length[node] = 0
        for node in range(size - 1, 0, -1):
            self._pull(node)

    def _set_leaf(self, node: int, value: int) -> None:
        self.length[node] = 1
        self.total[node] = value
        if value < 0:
            self.min_prefix[node], self.argmin[node] = value, 1
        else:
            self.min_prefix[node], self.argmin[node] = 0, 0
        if value > 0:
            self.max_prefix[node], self.argmax[node] = value, 1
        else:
            self.max_prefix[node], self.argmax[node] = 0, 0

    @staticmethod
    def _earlier_best(
        first_value: int,
        first_position: int,
        second_value: int,
        second_position: int,
        *,
        minimum: bool,
    ) -> tuple[int, int]:
        if minimum:
            if second_value < first_value:
                return second_value, second_position
        elif second_value > first_value:
            return second_value, second_position
        if second_value == first_value and second_position < first_position:
            return second_value, second_position
        return first_value, first_position

    def _pull(self, node: int) -> None:
        left = 2 * node
        right = left + 1
        left_len = self.length[left]
        self.length[node] = left_len + self.length[right]
        self.total[node] = self.total[left] + self.total[right]

        right_min = self.total[left] + self.min_prefix[right]
        self.min_prefix[node], self.argmin[node] = self._earlier_best(
            self.min_prefix[left],
            self.argmin[left],
            right_min,
            left_len + self.argmin[right],
            minimum=True,
        )
        right_max = self.total[left] + self.max_prefix[right]
        self.max_prefix[node], self.argmax[node] = self._earlier_best(
            self.max_prefix[left],
            self.argmax[left],
            right_max,
            left_len + self.argmax[right],
            minimum=False,
        )

    def update(self, index: int, value: int) -> None:
        node = self.size + index
        self._set_leaf(node, value)
        node //= 2
        while node:
            self._pull(node)
            node //= 2

    def summary(self) -> tuple[int, int, int, int, int]:
        return (
            self.total[1],
            self.min_prefix[1],
            self.argmin[1],
            self.max_prefix[1],
            self.argmax[1],
        )


def _add_event(events: list[list[tuple[int, bool]]], h: int, column: int, state: bool) -> None:
    if h:
        events[h].append((column, state))


def optimized_extrema(q: int, m: int, c: int) -> Extrema:
    """Compute exact extrema in O(q log q) time and O(q) memory.

    Positions are divided into q blocks.  Blocks are generated in increasing
    phase ``h = a*m_star (mod q)`` so each column has only two membership
    events.  Their summaries are stored by spatial block ``a`` and then
    concatenated in increasing spatial order.  This avoids assuming anything
    about the BBL prediction or stopping at a width threshold.
    """

    p, m_star, _ell, t = knot_parameters(q, m, c)
    m %= p
    s0 = m_star % q
    if not s0:
        raise AssertionError("m_star is not a unit modulo q")

    states = [False] * q
    values = [0] * q
    events: list[list[tuple[int, bool]]] = [[] for _ in range(q)]

    # Columns b=1,...,q-1.  For x=aq+b, membership is
    # ((floor(m_star*b/q) + h) mod q) < t.
    for b in range(1, q):
        column = b - 1
        shift = (m_star * b // q) % q
        states[column] = shift < t
        values[column] = t - q if states[column] else t
        start = (-shift) % q
        end = (t - shift) % q
        _add_event(events, start, column, True)
        _add_event(events, end, column, False)

    # Column b=q consists of the block endpoints.  The least-positive
    # convention makes membership 1 <= (h+s0 mod q) <= t.
    column = q - 1
    states[column] = 1 <= s0 <= t
    values[column] = t - q if states[column] else t
    start = (1 - s0) % q
    end = (t + 1 - s0) % q
    _add_event(events, start, column, True)
    _add_event(events, end, column, False)

    tree = _SummaryTree(values)
    blocks: list[tuple[int, int, int, int, int] | None] = [None] * q
    m0 = m % q
    for h in range(q):
        if h:
            for event_column, new_state in events[h]:
                states[event_column] = new_state
                tree.update(event_column, t - q if new_state else t)
        a = (h * m0) % q
        if blocks[a] is not None:
            raise AssertionError("phase-to-block map is not a permutation")
        blocks[a] = tree.summary()

    offset = 0
    min_g = max_g = 0
    argmin = argmax = 0
    for a, optional_block in enumerate(blocks):
        if optional_block is None:
            raise AssertionError("missing block summary")
        total, block_min, block_argmin, block_max, block_argmax = optional_block
        min_value = offset + block_min
        min_position = a * q + block_argmin
        if min_position < p and (min_value < min_g or (min_value == min_g and min_position < argmin)):
            min_g, argmin = min_value, min_position
        max_value = offset + block_max
        max_position = a * q + block_argmax
        if max_position < p and (max_value > max_g or (max_value == max_g and max_position < argmax)):
            max_g, argmax = max_value, max_position
        offset += total
    if offset != 0:
        raise AssertionError("optimized walk has nonzero total increment")
    return Extrema(min_g, max_g, argmin, argmax)


def surface_genus(width: int, order: int) -> int:
    """Convert jump width to the genus of the taut order-``order`` surface.

    GL22, Proposition 5.1 and Theorem 1.8 give
    ``width(g) = (order-1) + breadth(Delta)`` and
    ``breadth(Delta) = 1-chi(F)``.  With ``order`` boundary components,
    ``chi(F)=2-2*genus-order``.
    """

    numerator = width - 2 * order + 2
    if numerator < 0 or numerator % 2:
        raise ValueError(
            f"width {width} does not give a nonnegative integral genus for order {order}"
        )
    return numerator // 2
