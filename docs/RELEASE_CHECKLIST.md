# Local release checklist

- [x] Literal forward-multiple reference walk with (g(0)=0).
- [x] Independent inverse-membership event-sweep implementation.
- [x] (L(25,6)) example reproduced.
- [x] All unit presentations through (q=15) cross-checked.
- [x] Family sweep through (q=30) cross-checked.
- [x] Divisor generators compared with direct equation scans.
- [x] Core-exchange and mirror residue maps recorded.
- [x] Native implementation cross-checked against direct Python.
- [x] Inclusive (2 <= q <= 600) sweep complete.
- [x] Every expected key and stored column verified.
- [x] Both directional set differences reported.
- [x] Manifest and SHA-256 checksums generated.

The final four boxes were checked after `g0sv verify --q-max 600` returned
`status=complete` with zero missing keys, extra keys, column errors, or set
differences.
