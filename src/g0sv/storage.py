"""Atomic partitioned-Parquet sweep and release verification."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import asdict
from datetime import datetime, timezone
from math import gcd
from pathlib import Path
from typing import Any, Iterable

from .core import optimized_extrema, surface_genus
from .families import (
    FamilyRecord,
    family_records,
    independent_family_records,
    independent_prediction_records,
    independent_presentation_predictions,
    prediction_records,
    presentation_predictions,
)

SCHEMA_VERSION = 1
COMPRESSION = "zstd"
COMPRESSION_LEVEL = 12
ROW_GROUP_SIZE = 65_536
ROOT = Path(__file__).resolve().parents[2]
NATIVE_SOURCE = ROOT / "native" / "verifier.cpp"
COMPUTE_SOURCES = (
    ROOT / "src" / "g0sv" / "core.py",
    ROOT / "src" / "g0sv" / "families.py",
    NATIVE_SOURCE,
)
ROW_DTYPE_FIELDS = (
    ("q", "<i4"),
    ("m", "<i4"),
    ("c", "<i4"),
    ("min_g", "<i4"),
    ("max_g", "<i4"),
    ("argmin", "<i4"),
    ("argmax", "<i4"),
    ("width", "<i4"),
    ("genus", "<i4"),
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_bytes(_json_bytes(value))
    os.replace(temporary, path)


def source_identity() -> dict[str, Any]:
    source_hashes = {str(path.relative_to(ROOT)): sha256_file(path) for path in COMPUTE_SOURCES}
    configuration = {
        "schema_version": SCHEMA_VERSION,
        "normalization": "least-positive-residue; g(0)=0; extrema over 0<=x<q^2",
        "compression": COMPRESSION,
        "compression_level": COMPRESSION_LEVEL,
        "row_group_size": ROW_GROUP_SIZE,
        "q_min": 2,
        "source_hashes": source_hashes,
    }
    configuration_id = hashlib.sha256(_json_bytes(configuration)).hexdigest()
    return {**configuration, "configuration_id": configuration_id}


def _imports():
    try:
        import numpy as np
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as error:  # pragma: no cover - exercised by installation failures
        raise RuntimeError(
            "The sweep requires numpy and pyarrow. Install with `python -m pip install -e '.[test]'`."
        ) from error
    return np, pa, pq


def build_native(build_directory: Path | None = None) -> Path:
    identity = source_identity()["source_hashes"][str(NATIVE_SOURCE.relative_to(ROOT))][:16]
    build_directory = build_directory or ROOT / "build"
    build_directory.mkdir(parents=True, exist_ok=True)
    executable = build_directory / f"g0sv-native-{identity}"
    if executable.exists():
        return executable
    compiler = os.environ.get("CXX") or shutil.which("clang++") or shutil.which("g++")
    if not compiler:
        raise RuntimeError("No C++17 compiler was found (tried CXX, clang++, and g++)")
    temporary = executable.with_name(f".{executable.name}.{os.getpid()}.tmp")
    command = [
        compiler,
        "-O3",
        "-DNDEBUG",
        "-std=c++17",
        "-pthread",
        str(NATIVE_SOURCE),
        "-o",
        str(temporary),
    ]
    subprocess.run(command, cwd=ROOT, check=True)
    temporary.chmod(0o755)
    os.replace(temporary, executable)
    return executable


def expected_key_count(q: int, presentations: Iterable[int] | None = None) -> int:
    presentation_count = len(list(presentations)) if presentations is not None else len(presentation_predictions(q))
    unit_count = sum(gcd(c, q) == 1 for c in range(1, q))
    return presentation_count * unit_count


def key_digest(q_values, m_values, c_values) -> str:
    digest = hashlib.sha256()
    for q, m, c in zip(q_values, m_values, c_values, strict=True):
        digest.update(int(q).to_bytes(4, "little", signed=True))
        digest.update(int(m).to_bytes(4, "little", signed=True))
        digest.update(int(c).to_bytes(4, "little", signed=True))
    return digest.hexdigest()


def expected_key_digest(q: int, presentations: Iterable[int] | None = None) -> str:
    m_values = sorted(presentations if presentations is not None else independent_presentation_predictions(q))
    units = [c for c in range(1, q) if gcd(c, q) == 1]
    return key_digest(
        (q for _m in m_values for _c in units),
        (m for m in m_values for _c in units),
        (c for _m in m_values for c in units),
    )


def _native_rows(q: int, presentations: list[int], threads: int, work_directory: Path):
    np, _pa, _pq = _imports()
    executable = build_native()
    work_directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"q{q}-", dir=work_directory) as temp_name:
        temp = Path(temp_name)
        input_path = temp / "presentations.txt"
        output_path = temp / "rows.bin"
        input_path.write_text("".join(f"{m}\n" for m in presentations))
        result = subprocess.run(
            [str(executable), str(q), str(input_path), str(output_path), str(threads), "1"],
            cwd=ROOT,
            check=True,
            text=True,
            capture_output=True,
        )
        rows = np.fromfile(output_path, dtype=np.dtype(list(ROW_DTYPE_FIELDS)))
    return rows, result.stderr.strip()


def _python_rows(q: int, presentations: list[int]):
    np, _pa, _pq = _imports()
    units = [c for c in range(1, q) if gcd(c, q) == 1]
    rows = np.empty(len(presentations) * len(units), dtype=np.dtype(list(ROW_DTYPE_FIELDS)))
    index = 0
    for m in presentations:
        for c in units:
            extrema = optimized_extrema(q, m, c)
            rows[index] = (
                q,
                m,
                c,
                extrema.min_g,
                extrema.max_g,
                extrema.argmin,
                extrema.argmax,
                extrema.width,
                surface_genus(extrema.width, q),
            )
            index += 1
    return rows, "python event-sweep engine"


def _partition_paths(output: Path, q: int) -> tuple[Path, Path]:
    directory = output / "extrema" / f"q={q}"
    return directory / "part-00000.parquet", directory / "partition.json"


def validate_partition(output: Path, q: int, *, full_keys: bool = True) -> dict[str, Any]:
    np, _pa, pq = _imports()
    parquet_path, metadata_path = _partition_paths(output, q)
    if not parquet_path.exists() or not metadata_path.exists():
        raise RuntimeError(f"q={q}: partition or metadata is missing")
    metadata = json.loads(metadata_path.read_text())
    checksum = sha256_file(parquet_path)
    if checksum != metadata.get("sha256"):
        raise RuntimeError(f"q={q}: partition checksum mismatch")
    parquet_metadata = pq.read_metadata(parquet_path)
    if parquet_metadata.num_rows != metadata.get("row_count"):
        raise RuntimeError(f"q={q}: Parquet row count disagrees with partition metadata")
    independent_presentations = independent_presentation_predictions(q)
    expected_rows = expected_key_count(q, independent_presentations)
    if parquet_metadata.num_rows != expected_rows:
        raise RuntimeError(f"q={q}: expected {expected_rows} rows, found {parquet_metadata.num_rows}")
    if full_keys:
        table = pq.read_table(parquet_path, columns=["q", "m", "c"])
        actual_digest = key_digest(
            table["q"].to_numpy(), table["m"].to_numpy(), table["c"].to_numpy()
        )
        expected_digest = expected_key_digest(q, independent_presentations)
        if actual_digest != expected_digest or actual_digest != metadata.get("key_sha256"):
            raise RuntimeError(f"q={q}: key coverage/order digest mismatch")
    return metadata


def write_partition(
    output: Path,
    q: int,
    *,
    engine: str,
    threads: int,
) -> dict[str, Any]:
    np, pa, pq = _imports()
    predictions = presentation_predictions(q)
    presentations = sorted(predictions)
    independent = independent_presentation_predictions(q)
    if predictions != independent:
        raise RuntimeError(f"q={q}: production family generators disagree with independent scan")
    if engine == "native":
        rows, engine_report = _native_rows(q, presentations, threads, output / "_work")
    elif engine == "python":
        rows, engine_report = _python_rows(q, presentations)
    else:
        raise ValueError(f"unknown engine: {engine}")

    expected_rows = expected_key_count(q, presentations)
    if len(rows) != expected_rows:
        raise RuntimeError(f"q={q}: engine emitted {len(rows)} rows, expected {expected_rows}")
    units = [c for c in range(1, q) if gcd(c, q) == 1]
    expected_m = np.repeat(np.asarray(presentations, dtype=np.int32), len(units))
    expected_c = np.tile(np.asarray(units, dtype=np.int32), len(presentations))
    if not np.array_equal(rows["m"], expected_m) or not np.array_equal(rows["c"], expected_c):
        raise RuntimeError(f"q={q}: native row order/key set mismatch")
    if not np.all(rows["q"] == q):
        raise RuntimeError(f"q={q}: engine emitted an incorrect q value")
    expected_is_bbl = np.fromiter(
        (int(c) in predictions[int(m)] for m, c in zip(rows["m"], rows["c"], strict=True)),
        dtype=np.bool_,
        count=len(rows),
    )

    table = pa.table(
        {
            "q": pa.array(rows["q"], type=pa.int16()),
            "p": pa.array(np.full(len(rows), q * q, dtype=np.int32), type=pa.int32()),
            "m": pa.array(rows["m"], type=pa.int32()),
            "c": pa.array(rows["c"], type=pa.int16()),
            "min_g": pa.array(rows["min_g"], type=pa.int32()),
            "max_g": pa.array(rows["max_g"], type=pa.int32()),
            "argmin": pa.array(rows["argmin"], type=pa.int32()),
            "argmax": pa.array(rows["argmax"], type=pa.int32()),
            "width": pa.array(rows["width"], type=pa.int32()),
            "genus": pa.array(rows["genus"], type=pa.int32()),
            "is_bbl": pa.array(expected_is_bbl, type=pa.bool_()),
        }
    )
    parquet_path, metadata_path = _partition_paths(output, q)
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = parquet_path.with_name(f".{parquet_path.name}.{os.getpid()}.tmp")
    pq.write_table(
        table,
        temporary,
        compression=COMPRESSION,
        compression_level=COMPRESSION_LEVEL,
        row_group_size=ROW_GROUP_SIZE,
        use_dictionary=["q", "p", "m", "c", "genus", "is_bbl"],
        write_statistics=True,
    )
    os.replace(temporary, parquet_path)
    metadata = {
        "q": q,
        "row_count": len(rows),
        "presentation_count": len(presentations),
        "unit_count": len(units),
        "sha256": sha256_file(parquet_path),
        "key_sha256": expected_key_digest(q, presentations),
        "configuration_id": source_identity()["configuration_id"],
        "engine": engine,
        "engine_report": engine_report,
        "completed_at": utc_now(),
    }
    atomic_json(metadata_path, metadata)
    validate_partition(output, q)
    return metadata


def write_presentation_index(output: Path, q_values: Iterable[int]) -> Path:
    _np, pa, pq = _imports()
    records = [record for q in q_values for record in prediction_records(q)]
    table = pa.table(
        {
            "q": pa.array([record.q for record in records], type=pa.int16()),
            "p": pa.array([record.q * record.q for record in records], type=pa.int32()),
            "m": pa.array([record.m for record in records], type=pa.int32()),
            "source_m": pa.array([record.source_m for record in records], type=pa.int32()),
            "theorem_case": pa.array([record.theorem_case for record in records], type=pa.int8()),
            "branch": pa.array([record.branch for record in records], type=pa.string()),
            "d": pa.array([record.d for record in records], type=pa.int16()),
            "epsilon": pa.array([record.epsilon for record in records], type=pa.int8()),
            "n": pa.array([record.n for record in records], type=pa.int16()),
            "coordinate_map": pa.array([record.coordinate_map for record in records], type=pa.string()),
            "c_multiplier": pa.array([record.c_multiplier for record in records], type=pa.int16()),
            "predicted_cs": pa.array([list(record.predicted_cs) for record in records], type=pa.list_(pa.int16())),
        }
    )
    path = output / "presentation_index.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    pq.write_table(table, temporary, compression=COMPRESSION, compression_level=COMPRESSION_LEVEL)
    os.replace(temporary, path)
    return path


def load_state(output: Path) -> dict[str, Any] | None:
    path = output / "run_state.json"
    return json.loads(path.read_text()) if path.exists() else None


def run_sweep(
    *,
    q_min: int,
    q_max: int,
    output: Path,
    engine: str,
    threads: int,
    resume: bool,
) -> dict[str, Any]:
    if q_min < 2 or q_max < q_min:
        raise ValueError("require 2 <= q_min <= q_max")
    identity = source_identity()
    state_path = output / "run_state.json"
    state = load_state(output)
    if state:
        if not resume:
            raise RuntimeError(f"{state_path} already exists; pass --resume")
        if state.get("configuration_id") != identity["configuration_id"]:
            raise RuntimeError(
                "The existing run was created by different compute sources/configuration; "
                "use a new output directory."
            )
    else:
        state = {
            "configuration_id": identity["configuration_id"],
            "configuration": identity,
            "created_at": utc_now(),
            "completed": {},
        }
    state["requested_q_min"] = q_min
    state["requested_q_max"] = q_max
    state["engine"] = engine
    state["threads"] = threads
    state["resume_command"] = (
        f"g0sv sweep --q-min {q_min} --q-max {q_max} --output {output} "
        f"--engine {engine} --threads {threads} --resume"
    )
    output.mkdir(parents=True, exist_ok=True)
    atomic_json(state_path, state)

    started = time.perf_counter()
    for q in range(q_min, q_max + 1):
        key = str(q)
        if key in state["completed"]:
            validate_partition(output, q)
            print(f"q={q}: validated existing partition", flush=True)
            continue
        partition_started = time.perf_counter()
        metadata = write_partition(output, q, engine=engine, threads=threads)
        metadata["wall_seconds"] = time.perf_counter() - partition_started
        state["completed"][key] = metadata
        state["updated_at"] = utc_now()
        atomic_json(state_path, state)
        print(
            f"q={q}: {metadata['row_count']} rows in {metadata['wall_seconds']:.3f}s",
            flush=True,
        )
    completed_q = sorted(int(q) for q in state["completed"])
    write_presentation_index(output, completed_q)
    state["sweep_wall_seconds_this_invocation"] = time.perf_counter() - started
    state["updated_at"] = utc_now()
    atomic_json(state_path, state)
    return state


def _git_info() -> dict[str, Any]:
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
        ).stdout.strip()
    except subprocess.CalledProcessError:
        revision = None
    status = subprocess.run(
        ["git", "status", "--short"], cwd=ROOT, text=True, capture_output=True, check=False
    ).stdout
    return {"revision": revision, "dirty": bool(status.strip())}


def _family_record_dict(record: FamilyRecord) -> dict[str, Any]:
    result = asdict(record)
    result["family"] = record.family
    result["predicted_cs"] = list(record.predicted_cs)
    return result


def verify_release(output: Path, q_min: int, q_max: int) -> dict[str, Any]:
    np, _pa, pq = _imports()
    reports_directory = output / "reports"
    reports_directory.mkdir(parents=True, exist_ok=True)
    coverage_missing: list[dict[str, int]] = []
    coverage_extra: list[dict[str, int]] = []
    computed_not_expected: list[dict[str, int]] = []
    expected_not_computed: list[dict[str, int]] = []
    column_errors: list[dict[str, Any]] = []
    per_q: list[dict[str, Any]] = []
    family_totals: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "provenance_rows": 0,
            "presentation_keys": set(),
            "tested_memberships": 0,
            "computed_genus_zero_memberships": 0,
            "predicted_memberships": 0,
        }
    )
    total_rows = total_presentations = total_computed = total_expected = total_intersection = 0

    for q in range(q_min, q_max + 1):
        validate_partition(output, q)
        parquet_path, _metadata_path = _partition_paths(output, q)
        table = pq.read_table(parquet_path)
        columns = {name: table[name].to_numpy(zero_copy_only=False) for name in table.column_names}
        independent_predictions = independent_presentation_predictions(q)
        expected_presentations = sorted(independent_predictions)
        expected_units = [c for c in range(1, q) if gcd(c, q) == 1]
        expected_keys = {(q, m, c) for m in expected_presentations for c in expected_units}
        actual_keys = {
            (int(row_q), int(m), int(c))
            for row_q, m, c in zip(columns["q"], columns["m"], columns["c"], strict=True)
        }
        coverage_missing.extend(
            {"q": key[0], "m": key[1], "c": key[2]} for key in sorted(expected_keys - actual_keys)
        )
        coverage_extra.extend(
            {"q": key[0], "m": key[1], "c": key[2]} for key in sorted(actual_keys - expected_keys)
        )

        width_expected = columns["max_g"].astype(np.int64) - columns["min_g"].astype(np.int64)
        genus_numerator = columns["width"].astype(np.int64) - 2 * q + 2
        bad_width = np.flatnonzero(columns["width"] != width_expected)
        bad_genus = np.flatnonzero(
            (genus_numerator < 0)
            | (genus_numerator % 2 != 0)
            | (columns["genus"].astype(np.int64) != genus_numerator // 2)
        )
        bad_positions = np.flatnonzero(
            (columns["argmin"] < 0)
            | (columns["argmin"] >= q * q)
            | (columns["argmax"] < 0)
            | (columns["argmax"] >= q * q)
        )
        for kind, indices in (
            ("width", bad_width),
            ("genus", bad_genus),
            ("position", bad_positions),
        ):
            for index in indices[:1000]:
                column_errors.append(
                    {
                        "kind": kind,
                        "q": q,
                        "m": int(columns["m"][index]),
                        "c": int(columns["c"][index]),
                    }
                )

        computed_by_m: dict[int, set[int]] = {m: set() for m in expected_presentations}
        stored_bbl_by_m: dict[int, set[int]] = {m: set() for m in expected_presentations}
        for m, c, width, is_bbl in zip(
            columns["m"], columns["c"], columns["width"], columns["is_bbl"], strict=True
        ):
            if int(width) == 2 * q - 2:
                computed_by_m[int(m)].add(int(c))
            if bool(is_bbl):
                stored_bbl_by_m[int(m)].add(int(c))

        q_computed = q_expected = q_intersection = 0
        for m in expected_presentations:
            computed = computed_by_m[m]
            expected = independent_predictions[m]
            stored = stored_bbl_by_m[m]
            if stored != expected:
                for c in sorted(stored ^ expected):
                    column_errors.append({"kind": "is_bbl", "q": q, "m": m, "c": c})
            q_computed += len(computed)
            q_expected += len(expected)
            q_intersection += len(computed & expected)
            computed_not_expected.extend(
                {"q": q, "m": m, "c": c} for c in sorted(computed - expected)
            )
            expected_not_computed.extend(
                {"q": q, "m": m, "c": c} for c in sorted(expected - computed)
            )

        for record in independent_prediction_records(q):
            summary = family_totals[record.family]
            summary["provenance_rows"] += 1
            summary["presentation_keys"].add((q, record.m))
            summary["tested_memberships"] += len(expected_units)
            summary["computed_genus_zero_memberships"] += len(computed_by_m[record.m])
            summary["predicted_memberships"] += len(record.predicted_cs)

        q_summary = {
            "q": q,
            "presentations": len(expected_presentations),
            "rows": len(table),
            "computed_genus_zero": q_computed,
            "expected_bbl": q_expected,
            "intersection": q_intersection,
            "computed_not_expected": q_computed - q_intersection,
            "expected_not_computed": q_expected - q_intersection,
        }
        per_q.append(q_summary)
        total_rows += len(table)
        total_presentations += len(expected_presentations)
        total_computed += q_computed
        total_expected += q_expected
        total_intersection += q_intersection

    family_summary: dict[str, Any] = {}
    for family, summary in sorted(family_totals.items()):
        presentation_keys = summary.pop("presentation_keys")
        family_summary[family] = {**summary, "distinct_presentations": len(presentation_keys)}

    complete = not (
        coverage_missing
        or coverage_extra
        or computed_not_expected
        or expected_not_computed
        or column_errors
    )
    summary = {
        "status": "complete" if complete else "incomplete",
        "q_min": q_min,
        "q_max": q_max,
        "inclusive": True,
        "deduplicated": {
            "presentations": total_presentations,
            "rows": total_rows,
            "computed_genus_zero": total_computed,
            "expected_bbl": total_expected,
            "intersection": total_intersection,
            "computed_not_expected": len(computed_not_expected),
            "expected_not_computed": len(expected_not_computed),
        },
        "coverage": {
            "missing_keys": len(coverage_missing),
            "extra_keys": len(coverage_extra),
            "column_errors": len(column_errors),
        },
        "per_q": per_q,
        "per_family_provenance_weighted": family_summary,
        "verified_at": utc_now(),
    }
    discrepancies = {
        "computed_not_expected": computed_not_expected,
        "expected_not_computed": expected_not_computed,
        "paper_readme_coordinate_note": (
            "The README's direct case-2 row for (q,m,d,epsilon,n)=(7,30,5,+1,3) predicts "
            "only c={1,2,5,6}. Paper-authoritative core exchange from the coincident case-3 "
            "n=2 presentation m=18 sends c to 30c=2c (mod 7), transporting its exceptional "
            "+/-2 pair to c=+/-3={3,4} at m=30. The certified expected set is therefore all "
            "six units. Separately, paper Theorem 1.4's signed phrase q=1-2d becomes q=2d+/-1 "
            "under Theorems 6.1 and 10.17, exactly README case 3 with n=2."
        ),
    }
    coverage = {
        "missing_keys": coverage_missing,
        "extra_keys": coverage_extra,
        "column_errors": column_errors,
    }
    atomic_json(reports_directory / "summary.json", summary)
    atomic_json(reports_directory / "coverage.json", coverage)
    atomic_json(reports_directory / "discrepancies.json", discrepancies)

    report_lines = [
        "# Verification summary",
        "",
        f"Status: **{summary['status']}**",
        "",
        f"Inclusive range: `{q_min} <= q <= {q_max}`",
        "",
        f"Deduplicated presentations: {total_presentations:,}",
        "",
        f"Per-knot rows: {total_rows:,}",
        "",
        f"Computed genus-zero memberships: {total_computed:,}",
        "",
        f"Expected BBL memberships: {total_expected:,}",
        "",
        f"Intersection: {total_intersection:,}",
        "",
        f"Computed but not expected: {len(computed_not_expected):,}",
        "",
        f"Expected but not computed: {len(expected_not_computed):,}",
        "",
        f"Missing keys: {len(coverage_missing):,}; extra keys: {len(coverage_extra):,}; column errors: {len(column_errors):,}.",
        "",
        "The README's direct case-2 row at `(q,m,d,epsilon,n)=(7,30,5,+1,3)` omits the transported exceptional pair `c=+/-3`. Paper-authoritative core exchange from the coincident case-3 `n=2` presentation `m=18` sends `c` to `30c=2c (mod 7)`, so the certified expected set at `m=30` is all six units.",
        "",
        "The signed-parameter phrase `q = 1 - 2d` in paper Theorem 1.4 becomes `q = 2d +/- 1` under the paper's positive-q normalization (Theorems 6.1 and 10.17), exactly case 3 with `n = 2`.",
        "",
    ]
    (reports_directory / "summary.md").write_text("\n".join(report_lines))

    index_path = write_presentation_index(output, range(q_min, q_max + 1))
    files = [
        path
        for path in output.rglob("*")
        if path.is_file()
        and path.name != "manifest.json"
        and "_work" not in path.parts
    ]
    checksums = {str(path.relative_to(output)): sha256_file(path) for path in sorted(files)}
    manifest = {
        "status": summary["status"],
        "dataset_license": "CC0-1.0",
        "code_license": "MIT",
        "q_min": q_min,
        "q_max": q_max,
        "inclusive": True,
        "row_count": total_rows,
        "presentation_count": total_presentations,
        "configuration": source_identity(),
        "git": _git_info(),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
        },
        "presentation_index": str(index_path.relative_to(output)),
        "checksums_sha256": checksums,
        "created_at": utc_now(),
    }
    atomic_json(output / "manifest.json", manifest)
    return {"summary": summary, "coverage": coverage, "discrepancies": discrepancies, "manifest": manifest}
