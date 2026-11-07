# Mathematical and coordinate specification

The authority for the BBL classes is *Genus Zero Simple Knots* by Joshua Evan
Greene, Shiyu Liang, and John Luecke (2026), especially Theorems 1.4, 6.1,
10.17, and 10.18 and Section 10.1. The implementation uses the paper's
least-positive-residue convention, $m^{\mathord{\ast}}$ with
$mm^{\mathord{\ast}}\equiv1\pmod{q^2}$, and the normalization $g(0)=0$.

For $k=cq$, $1\le c<q$, and $\gcd(c,q)=1$, put

$$
\ell=[m^{\mathord{\ast}}k]_{q^2},\qquad t=\ell/q,
\qquad S=\lbrace[jm]_{q^2}:1\le j\le\ell\rbrace.
$$

Then

$$
g(x)=g(x-1)+t-q\mathbf{1}_S(x),\qquad 1\le x\le q^2.
$$

Extrema and their attaining positions are taken on $0\le x<q^2$; the value at
$x=q^2$ is the duplicate $g(0)=0$. The implementation records the earliest
position if the same extremal value occurs more than once.

## Genus conversion

[Greene-Luecke, *Fibered simple knots*, Proposition 5.1 and Theorem
1.8](https://arxiv.org/abs/2106.08485) give

$$
\mathrm{width}(g)=(q-1)+\mathrm{breadth}(\Delta_K),
\qquad \mathrm{breadth}(\Delta_K)=1-\chi(F).
$$

The taut rational Seifert surface has $q$ boundary components, so
$\chi(F)=2-2\mathrm{genus}(F)-q$. Therefore the integer stored in the
`genus` column is

$$
\mathrm{genus}(F)=
\frac{\mathrm{width}(g)-2q+2}{2}.
$$

The verifier checks that this value is integral and nonnegative but does not
use it to compute the width.

## Precise discrepancy at q=7

The README's direct case-2 row for
$(q,m,d,\epsilon,n)=(7,30,5,+1,3)$ gives only
$c=\pm1,\pm5=\lbrace1,2,5,6\rbrace$. This is not the full paper-authoritative set in
that marked presentation.

The paper also gives the coincident case-3 $n=2$ presentation $L(49,18)$,
where $18\equiv30^{-1}\pmod{49}$. Its exceptional $c=\pm2$ pair must be
transported by the paper's core exchange. Since source $\mu$ maps to target
$\mu'=30\mu$, the residue map in the $L(49,30)$ marking is

$$
c\longmapsto30c\equiv2c\pmod7.
$$

Thus the exceptional pair becomes $c=\pm3=\lbrace3,4\rbrace$, and the certified
expected set at $L(49,30)$ is all six units modulo 7. The presentation index
retains `source_m`, `coordinate_map`, and `c_multiplier` so this addition is
auditable rather than silently folded into the prediction.

Paper Theorem 1.4 also phrases its exceptional condition in signed parameters
as $q=1-2d$. Under the normalization with positive $q$ in Theorems 6.1 and 10.17,
this is $q=2d\pm1$, exactly case 3 with $n=2$.
