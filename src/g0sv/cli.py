"""Command-line interface for the verifier."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from math import gcd
from pathlib import Path

from .core import optimized_extrema, reference_walk, reference_walk_points, surface_genus
from .families import presentation_predictions
from .storage import (
    ROOT,
    ROW_DTYPE_FIELDS,
    _imports,
    build_native,
    run_sweep,
    verify_release,
)


def _walk(args: argparse.Namespace) -> int:
    if args.stream:
        writer = csv.writer(sys.stdout)
        writer.writerow(("x", "g", "increment", "in_s"))
        for point in reference_walk_points(args.q, args.m, args.c):
            writer.writerow(point)
        return 0
    extrema = (
        reference_walk(args.q, args.m, args.c)
        if args.engine == "reference"
        else optimized_extrema(args.q, args.m, args.c)
    )
    p = args.q * args.q
    m = args.m % p
    m_star = pow(m, -1, p)
    ell = (m_star * args.c * args.q) % p or p
    result = {
        "q": args.q,
        "p": p,
        "m": m,
        "m_star": m_star,
        "c": args.c,
        "ell": ell,
        "t": ell // args.q,
        "min_g": extrema.min_g,
        "max_g": extrema.max_g,
        "argmin": extrema.argmin,
        "argmax": extrema.argmax,
        "width": extrema.width,
        "genus": surface_genus(extrema.width, args.q),
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


def _sweep(args: argparse.Namespace) -> int:
    output = args.output.resolve()
    try:
        state = run_sweep(
            q_min=args.q_min,
            q_max=args.q_max,
            output=output,
            engine=args.engine,
            threads=args.threads,
            resume=args.resume,
        )
    except KeyboardInterrupt:
        print("Sweep interrupted; completed partitions are intact.", file=sys.stderr)
        print(
            f"Resume with: g0sv sweep --q-min {args.q_min} --q-max {args.q_max} "
            f"--output {output} --engine {args.engine} --threads {args.threads} --resume",
            file=sys.stderr,
        )
        return 130
    print(json.dumps({"completed_q": sorted(map(int, state["completed"])), "output": str(output)}, indent=2))
    return 0


def _verify(args: argparse.Namespace) -> int:
    result = verify_release(args.output.resolve(), args.q_min, args.q_max)
    print(json.dumps(result["summary"]["deduplicated"], indent=2, sort_keys=True))
    print(f"status={result['summary']['status']}")
    return 0 if result["summary"]["status"] == "complete" else 1


def _benchmark_native(q: int, threads: int, repeat: int) -> dict:
    np, _pa, _pq = _imports()
    presentations = sorted(presentation_predictions(q))
    executable = build_native()
    with tempfile.TemporaryDirectory(prefix="g0sv-bench-") as name:
        directory = Path(name)
        inputs = directory / "presentations.txt"
        output = directory / "rows.bin"
        inputs.write_text("".join(f"{m}\n" for m in presentations))
        result = subprocess.run(
            [str(executable), str(q), str(inputs), str(output), str(threads), str(repeat)],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        rows = np.fromfile(output, dtype=np.dtype(list(ROW_DTYPE_FIELDS)))
    match = re.search(r"seconds=([0-9.eE+-]+) rows_per_second=([0-9.eE+-]+)", result.stderr)
    if not match:
        raise RuntimeError(f"could not parse native benchmark: {result.stderr}")
    return {
        "q": q,
        "presentations": len(presentations),
        "rows": len(rows),
        "repeat": repeat,
        "seconds": float(match.group(1)),
        "rows_per_second": float(match.group(2)),
    }


def _benchmark(args: argparse.Namespace) -> int:
    direct_cases = [(5, 6, 1), (15, 16, 7), (30, 31, 29)]
    direct_results = []
    for q, m, c in direct_cases:
        started = time.perf_counter()
        direct = reference_walk(q, m, c)
        direct_seconds = time.perf_counter() - started
        started = time.perf_counter()
        optimized = optimized_extrema(q, m, c)
        optimized_seconds = time.perf_counter() - started
        if direct != optimized:
            raise RuntimeError(f"benchmark cross-check failed for {(q, m, c)}")
        direct_results.append(
            {
                "q": q,
                "m": m,
                "c": c,
                "direct_seconds": direct_seconds,
                "optimized_python_seconds": optimized_seconds,
            }
        )
    native = [_benchmark_native(q, args.threads, args.repeat) for q in args.q_values]
    result = {
        "threads": args.threads,
        "direct_cross_checks": direct_results,
        "native_partitions": native,
        "note": "Full-sweep duration is measured by the sweep itself; throughput depends strongly on q.",
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


def _inspect(args: argparse.Namespace) -> int:
    _np, _pa, pq = _imports()
    path = args.output.resolve() / "extrema" / f"q={args.q}" / "part-00000.parquet"
    table = pq.read_table(path)
    rows = table.to_pylist()
    matches = [row for row in rows if row["m"] == args.m % (args.q * args.q) and row["c"] == args.c]
    if len(matches) != 1:
        raise RuntimeError(f"expected one stored row, found {len(matches)}")
    stored = matches[0]
    regenerated = optimized_extrema(args.q, args.m, args.c)
    result = {
        "stored": stored,
        "regenerated": {
            "min_g": regenerated.min_g,
            "max_g": regenerated.max_g,
            "argmin": regenerated.argmin,
            "argmax": regenerated.argmax,
            "width": regenerated.width,
            "genus": surface_genus(regenerated.width, args.q),
        },
    }
    result["matches"] = all(
        stored[key] == result["regenerated"][key]
        for key in ("min_g", "max_g", "argmin", "argmax", "width", "genus")
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["matches"] else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="g0sv")
    subparsers = parser.add_subparsers(dest="command", required=True)

    walk = subparsers.add_parser("walk", help="compute or stream one complete walk")
    walk.add_argument("--q", type=int, required=True)
    walk.add_argument("--m", type=int, required=True)
    walk.add_argument("--c", type=int, required=True)
    walk.add_argument("--engine", choices=("reference", "optimized"), default="reference")
    walk.add_argument("--stream", action="store_true", help="stream the literal walk as CSV")
    walk.set_defaults(func=_walk)

    sweep = subparsers.add_parser("sweep", help="write atomic resumable q partitions")
    sweep.add_argument("--q-min", type=int, default=2)
    sweep.add_argument("--q-max", type=int, default=600)
    sweep.add_argument("--output", type=Path, default=Path("release"))
    sweep.add_argument("--engine", choices=("native", "python"), default="native")
    sweep.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1))
    sweep.add_argument("--resume", action="store_true")
    sweep.set_defaults(func=_sweep)

    verify = subparsers.add_parser("verify", help="check coverage, columns, and both set differences")
    verify.add_argument("--q-min", type=int, default=2)
    verify.add_argument("--q-max", type=int, default=600)
    verify.add_argument("--output", type=Path, default=Path("release"))
    verify.set_defaults(func=_verify)

    benchmark = subparsers.add_parser("benchmark", help="cross-check and benchmark representative inputs")
    benchmark.add_argument("--q-values", type=int, nargs="+", default=[100, 300, 600])
    benchmark.add_argument("--repeat", type=int, default=3)
    benchmark.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1))
    benchmark.add_argument("--output", type=Path, default=Path("reports/benchmark.json"))
    benchmark.set_defaults(func=_benchmark)

    inspect = subparsers.add_parser("inspect", help="inspect a stored row and regenerate its extrema")
    inspect.add_argument("--output", type=Path, default=Path("release"))
    inspect.add_argument("--q", type=int, required=True)
    inspect.add_argument("--m", type=int, required=True)
    inspect.add_argument("--c", type=int, required=True)
    inspect.set_defaults(func=_inspect)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
