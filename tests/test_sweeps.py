from __future__ import annotations

import subprocess
from math import gcd

import numpy as np

from g0sv.core import optimized_extrema, reference_walk, surface_genus
from g0sv.families import presentation_predictions
from g0sv.storage import ROOT, ROW_DTYPE_FIELDS, build_native


def test_family_sweep_through_q30_against_reference() -> None:
    for q in range(2, 31):
        for m in presentation_predictions(q):
            for c in range(1, q):
                if gcd(c, q) == 1:
                    assert optimized_extrema(q, m, c) == reference_walk(q, m, c)


def test_native_engine_against_direct_walk(tmp_path) -> None:
    executable = build_native(tmp_path / "build")
    for q in (30, 100):
        presentations = sorted(presentation_predictions(q))
        input_path = tmp_path / f"m-{q}.txt"
        output_path = tmp_path / f"rows-{q}.bin"
        input_path.write_text("".join(f"{m}\n" for m in presentations))
        subprocess.run(
            [str(executable), str(q), str(input_path), str(output_path), "4", "1"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        rows = np.fromfile(output_path, dtype=np.dtype(list(ROW_DTYPE_FIELDS)))
        # All q=30 rows and a deterministic spread of larger q=100 rows.
        selected = rows if q == 30 else rows[:: max(1, len(rows) // 12)]
        for row in selected:
            expected = reference_walk(q, int(row["m"]), int(row["c"]))
            assert (
                int(row["min_g"]),
                int(row["max_g"]),
                int(row["argmin"]),
                int(row["argmax"]),
                int(row["width"]),
                int(row["genus"]),
            ) == (
                expected.min_g,
                expected.max_g,
                expected.argmin,
                expected.argmax,
                expected.width,
                surface_genus(expected.width, q),
            )
