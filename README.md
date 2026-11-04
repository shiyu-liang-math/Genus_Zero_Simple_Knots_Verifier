# Genus Zero Simple Knots · computational verification

Companion project to *Genus Zero Simple Knots* by Joshua Evan Greene, Shiyu
Liang, and John Luecke (2026). The paper classifies genus-zero simple knots in
lens spaces. This project tests the classification in the four families below
for **every integer $2\le q\le600$**, including composite $q$. The largest
lens-space order in the search is $q^2=360{,}000$.

**Verification status: complete for the inclusive range $2\le q\le600$.**
The release verifier checked 36,548,348 per-knot rows across 141,834 distinct
marked presentations. It found 510,740 computed genus-zero memberships and
510,740 paper-authoritative BBL memberships, with intersection 510,740. Both
directional set differences, both coverage differences, and the derived-column
error set are empty. See [the compact report](release/reports/summary.md).

## The comparison

Fix an ordered presentation $L(q^2,m)$ and the homology generator $\mu$
used in the paper. Write the knot class as $cq\mu$. Taking all
$1\le c<q$ with $\gcd(c,q)=1$ enumerates exactly the classes of order
$q$. Both signs are retained as individual residues.

For each marked presentation the verifier compares two independently generated
sets modulo $q$:

- **Computed:** every $c$ whose exact jump-function width is $2q-2$.
- **Expected:** the paper's BBL classes transported into the same marked
  coordinates, retaining mirror/core-exchange maps and unioning coincident
  descriptions.

The extrema engine never consults the expected set and never stops at a width
threshold.

### Paper/README coordinate correction at q=7

The raw case-2 table row for
$(q,m,d,\epsilon,n)=(7,30,5,+1,3)$ supplies
$c=\pm1,\pm5=\{1,2,5,6\}$, but this is not the full paper-authoritative
union in that marking. The coincident case-3 $n=2$ presentation
$L(49,18)$, with $18\equiv30^{-1}\pmod{49}$, contributes the exceptional
$c=\pm2$ pair. Under core exchange the marked residue map is
$c\mapsto30c\equiv2c\pmod7$, so that pair becomes
$c=\pm3=\{3,4\}$ at $m=30$. The certified expected set there is all six
units. See [the mathematical specification](MATHEMATICAL_SPEC.md) for the
full audit trail.

## Lens-space families

The positive-presentation generators are:

| Paper case and hypotheses | $m$ | BBL $c$ |
|---|:---:|---|
| Case 1: $\gcd(d,q)=1$<br>$1\le d\le q$ | $dq+1$ | $\pm1,\pm d$ |
| Case 1: $\gcd(d,q)=2$<br>$1\le d\le q$ | $dq+1$ | $\pm1$ |
| Case 2: $2q+\epsilon=nd$<br>$d,n\ge3$ odd | $d(q-\epsilon)$ | $\pm1,\pm d$ |
| Case 3: $q-\epsilon=nd$<br>$d\ge3$ odd, $n>1$ | $d(q-\epsilon)$ | $\pm1,\pm d$<br>$\pm2$ if $n=2$ |

Both signs $\epsilon=\pm1$ are used. The presentation index preserves every
family row and every induced coordinate map; the extrema table is deduplicated
on $(q,m\bmod q^2,c)$.

## Computing the extrema

Put $p=q^2$, $k=cq$, and choose $m^*$ with
$mm^*\equiv1\pmod p$. With $[x]_p$ denoting the least positive residue,
so $[0]_p=p$, define

$$
\ell=[m^*k]_p,\qquad t=\ell/q,\qquad
S=\{[jm]_p:1\le j\le\ell\}.
$$

Starting at $g(0)=0$, take

$$
g(x)=g(x-1)+t-q\mathbf1_S(x),\qquad1\le x\le p.
$$

The per-knot extrema table stores the exact minimum and maximum over
$0\le x<p$, together with attaining positions `argmin` and `argmax`. The
implementation selects the earliest position if an extremum is attained more
than once. The integer surface genus is independently derived as

$$
\mathrm{genus}(F)=\frac{\max g-\min g-2q+2}{2}.
$$

Its derivation and source are recorded in
[MATHEMATICAL_SPEC.md](MATHEMATICAL_SPEC.md).

## Install and test

Python 3.11+ and a C++17 compiler are required for the release sweep.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m pytest
```

The tests exhaust every unit $m\bmod q^2$ and admissible $c$ through
$q=15$, cross-check the complete family sweep through $q=30$, compare the
native worker with the literal Python walk, scan family equations independently,
and exercise the $q=2$ and $q=7$ boundaries.

On the release machine (8 worker threads), the native engine measured roughly
403k rows/s at $q=100$, 153k rows/s at $q=300$, and 76k rows/s at
$q=600$. The complete 36.55-million-row sweep, including atomic Parquet
writes and per-partition validation, took 455.25 seconds.

## Reproduce L(25,6)

```bash
g0sv walk --q 5 --m 6 --c 1 --engine reference
```

This reports $m^*=21$, $\ell=5$, $t=1$, unique maximum $g(4)=4$, unique
minimum $g(6)=-4$, and width $8=2q-2$.
Stream the complete normalized walk as CSV with `--stream`.

## Run, resume, and verify the release

```bash
g0sv benchmark --q-values 100 300 600 --threads 8
g0sv sweep --q-min 2 --q-max 600 --output release --engine native --threads 8
g0sv verify --q-min 2 --q-max 600 --output release
```

Partition writes are atomic. A stopped run is resumed without duplicating rows:

```bash
g0sv sweep --q-min 2 --q-max 600 --output release --engine native --threads 8 --resume
```

Resume validates the source/configuration identity, SHA-256 checksum, Parquet
row count, and independently generated key digest of every completed partition.

Inspect a stored row and regenerate it:

```bash
g0sv inspect --output release --q 5 --m 6 --c 1
```

## Release layout

- `release/extrema/q=<q>/part-00000.parquet`: Zstandard-compressed per-knot
  extrema, partitioned by $q$.
- `release/presentation_index.parquet`: family parameters and coordinate maps.
- `release/reports/coverage.json`: missing/extra keys and column validation.
- `release/reports/discrepancies.json`: both directional set differences with
  concrete $(q,m,c)$ records.
- `release/reports/summary.{json,md}`: per-$q$, per-family, and deduplicated
  totals.
- `release/manifest.json`: range, row counts, environment, code identity, Git
  revision, and SHA-256 checksums.

The verifier source is licensed under the [MIT License](LICENSE). Generated
data and machine-generated reports are dedicated under
[CC0 1.0 Universal](DATA_LICENSE.md).

## Citation request

If you use this verifier or its generated dataset, please cite
*Genus Zero Simple Knots* by **Joshua Evan Greene, Shiyu Liang, and John Luecke
(2026)**. Please also cite the software release after the author publishes it.
