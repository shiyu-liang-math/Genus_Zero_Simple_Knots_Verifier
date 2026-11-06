# Using the verification dataset

The release is a **partitioned Parquet dataset**, not a database server or a
SQLite file. Parquet is a compressed, typed, column-oriented format commonly
used for large analytical datasets. It is efficient for 36.5 million rows, but
it requires a Parquet reader such as PyArrow, DuckDB, pandas, R, Julia, or
Mathematica.

After extracting the release archive, the relevant layout is:

```text
release/
  extrema/
    q=2/
      part-00000.parquet   # the actual knot records
      partition.json       # integrity and resume metadata only
    q=3/
      part-00000.parquet
      partition.json
    ...
    q=600/
      part-00000.parquet
      partition.json
  presentation_index.parquet
  reports/
  manifest.json
  run_state.json
```

The `partition.json` files do not contain extrema. They record the row count,
checksums, engine, configuration identity, and completion information for the
corresponding `part-00000.parquet` file.

## Extrema columns

Each row of an extrema Parquet file represents one order-$q$ knot class in one
marked presentation. The schema is:

```text
q        int16   order of the knot class
p        int32   lens-space order, equal to q^2
m        int32   marked lens-space presentation L(p,m)
c        int16   unit residue with knot class k=cq
min_g    int32   exact minimum of g on 0 <= x < p
max_g    int32   exact maximum of g on 0 <= x < p
argmin   int32   earliest x attaining min_g
argmax   int32   earliest x attaining max_g
width    int32   max_g-min_g
genus    int32   (width-2q+2)/2
is_bbl   bool    true when c belongs to the predicted BBL set
```

The table stores one attaining position for each extremum. If the same value is
attained more than once, `argmin` or `argmax` is the earliest such position.

## Quickest Python example

Install PyArrow:

```bash
python3 -m pip install pyarrow
```

Print every record in the $q=2$ partition:

```bash
python3 - <<'PY'
import pyarrow.parquet as pq

table = pq.read_table("release/extrema/q=2/part-00000.parquet")
for row in table.to_pylist():
    print(row)
PY
```

Read one knot, for example $(q,m,c)=(5,6,1)$, without loading the full
dataset:

```bash
python3 - <<'PY'
import pyarrow.compute as pc
import pyarrow.parquet as pq

table = pq.read_table("release/extrema/q=5/part-00000.parquet")
mask = pc.and_(pc.equal(table["m"], 6), pc.equal(table["c"], 1))
print(table.filter(mask).to_pylist())
PY
```

The result contains
`min_g=-4`, `max_g=4`, `argmin=6`, `argmax=4`, `width=8`, and `genus=0`.

## Query with DuckDB

DuckDB provides SQL directly over the Parquet files. On macOS:

```bash
brew install duckdb
```

Query one knot across all partitions:

```bash
duckdb -c "
SELECT q, m, c, min_g, max_g, argmin, argmax, width, genus, is_bbl
FROM read_parquet(
  'release/extrema/q=*/part-00000.parquet',
  hive_partitioning=false
)
WHERE q=5 AND m=6 AND c=1;
"
```

List every computed genus-zero knot for $q=25$:

```bash
duckdb -c "
SELECT m, c, min_g, max_g, argmin, argmax
FROM read_parquet(
  'release/extrema/q=25/part-00000.parquet',
  hive_partitioning=false
)
WHERE genus=0
ORDER BY m, c;
"
```

Export a partition to an ordinary CSV file:

```bash
duckdb -c "
COPY (
  SELECT *
  FROM read_parquet(
    'release/extrema/q=5/part-00000.parquet',
    hive_partitioning=false
  )
  ORDER BY m, c
) TO 'q5-extrema.csv' (HEADER, DELIMITER ',');
"
```

## Inspect and independently regenerate one row

After installing the verifier from the repository, its `inspect` command reads
one stored row and independently recomputes its extrema:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/g0sv inspect --output release --q 5 --m 6 --c 1
```

The returned `matches` field is `true` when the regenerated and stored extrema
agree.

## Other release files

- `presentation_index.parquet` records every family description and coordinate
  transport used to construct the expected sets.
- `reports/coverage.json` reports missing keys, extra keys, and column errors.
- `reports/discrepancies.json` contains both directional set differences.
- `reports/summary.json` and `summary.md` contain aggregate totals.
- `manifest.json` records the range, schema configuration, licenses, source
  identity, environment, and SHA-256 checksums.
- `run_state.json` records the completed partition state used by resumable
  sweeps.
