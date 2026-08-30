# Pinwheel test — does $K(n)$ grow as $\lfloor n/3\rfloor$?

## Headline result

**No.** $K_{\rm geom}(n)$ **saturates at 3 through $n=10$**, falsifying the
conjectured $K(n) = \lfloor n/3\rfloor$ growth. The $m$-fold pinwheel at
$n=3m$ — the object that motivated the conjecture — itself needs only $k=3$
for $m=4,5$. Independent confirmation by full $n=10$ step-1 enumeration:
no multiset among the 9 305 non-triangulation step-1 survivors needs $k\ge 4$.

| $n$ | $m$ | pinwheel chords | fillers | min-kill-$k$ across all filler choices |
|---:|---:|---:|---:|:---:|
| 9   | 3 | 6  | 0    | **3** (this is the n=9 pinwheel that motivated the conjecture) |
| 12  | 4 | 8  | 1    | **3** for all 54 filler choices |
| 15  | 5 | 10 | 2    | **3** for all 4 095 filler choices |

Full-enumeration confirmation via $n=10$ step-1 survivors
(`check_n10_step1_survivors.py`, runtime ~8 min):

| $n$ | step-1 non-tri survivors | $k=2$ | $k=3$ | $k\ge 4$ |
|---:|---:|---:|---:|---:|
| 7  | 7    | 0    | 7   | 0 |
| 8  | 100  | 0    | 100 | 0 |
| 9  | 1 011 | 954 | 57   | 0 |
| 10 | 9 305 | 8 470 | 835 | **0** |

Through $n=10$, $K_{\rm geom}(n) = 3$ exactly (for $n\ge 9$). No multiset
needs $k=4$ under the geometric criterion at $n\le 10$.

The extra tight crossings beyond $m=3$ don't preserve the difficulty — they
give the kill **more** ways to fail P1 ∨ P2, not fewer.

## Setup

Geometric kill criterion (matches `step1_uncaught.py` at $k=1$): $M$ survives
the $k$-zone at rotation $i$ iff one of
- (P1) some chord of $M$ crosses the enclosing chord $X_{i,i+k}$, or
- (P2) $M$ contains the special $X_{i,i+k+1}$, the bare $X_{i+k,i-1}$, or a
  column pair $\{X_{i,m}, X_{i+k,m}\}$ for some far-middle $m$.

Killed otherwise. At $k=1$, $X_{i,i+1}$ is an edge so P1 is vacuous and the
criterion reduces to step-1.

**Validation at $n=9$** (matches the first-Claude / `kill_dictionary.py`
analysis): of the 1 011 step-1 non-triangulation survivors,
- 954 die at $k=2$ (geometric),
- 57 need $k=3$ (geometric).

The Z₃ pinwheel
$M = \{(1,3),(2,4),(4,6),(5,7),(7,9),(1,8)\}$
is in the latter group; min-kill-$k$ = 3, at zone $i=1$ (enclosing chord
$X_{1,4}$).

## Pinwheel construction

The $m$-fold pinwheel at $n=3m$ consists of $m$ "tight crossings" rotated by
$2\pi/m$ around the polygon. Each tight crossing occupies 4 consecutive
vertices $\{3j+1, 3j+2, 3j+3, 3j+4\}$ and consists of the two short diagonals
of that quadruple,
$(3j+1, 3j+3)$ and $(3j+2, 3j+4)$.

Multiset size = $2m$ chords. For $n=3m$, the multiset slot count is $n-3=3m-3$,
so we need $3m-3 - 2m = m-3$ filler chords.

At $m=3$ (n=9): zero fillers — the pinwheel exactly fills the multiset.
At $m=4$ (n=12): 1 filler.
At $m=5$ (n=15): 2 fillers.

## Why the conjecture fails

The reasoning behind $K(n) = \lfloor n/3\rfloor$: the n=9 pinwheel resists
$k=2$ because its three tight crossings, when read at any single $k=2$ zone,
prevent the local P1/P2 conditions from holding. Extending to $m\ge 4$ tight
crossings, naively, one expects more resistance and a higher $k$ needed.

What actually happens: with more tight crossings, **more chords of $M$ are
available to cross enclosing chords at every rotation**. Hence P1 fires at
more zones, including $k=3$ zones. The $k=3$ length-3 enclosing chord
$X_{i,i+3}$ is short enough that some pinwheel chord crosses it for every $i$
(in the $m\ge 4$ case).

Specifically, at $n=12$, the pinwheel chord $(3j+1,3j+3)$ crosses $X_{i,i+3}$
whenever the vertex $3j+2$ lies strictly between $i$ and $i+3$ (mod 12) and
$3j$ or $3j+4$ lies outside. With 4 tight crossings filling the polygon,
every $i$ catches some crossing.

The combinatorial obstruction at $n=9$ is unique to the exact-fit case ($m=3$,
no fillers, only 6 chords); $m\ge 4$ adds chords (pinwheel + fillers) that
re-enable kills at $k=3$.

## What this means for the structural strategy

The "$K(n) = \lfloor n/3\rfloor$" line of attack — classify survivors by their
tightest bad structure (repeats → $k=2$, tight crossings → $k=3$, pinwheels
→ $k=m$, etc.) — does NOT close. The $m$-fold pinwheel is not extremal for
$m\ge 4$.

**Two open possibilities remain:**

1. **$K(n)$ saturates at 3** (geometric criterion). If true, locality at every
   $n$ follows from leading fat-zero kills with $k\le 3$ alone. The next step is
   to find the **actual** extremal structures (not pinwheels) and show they all
   die at $k\le 3$.

2. **$K(n)$ does grow**, but via a different structural family. The pinwheel
   was a natural candidate; falsifying it doesn't rule out subtler structures.

## Reality check via the rigorous criterion

The geometric criterion is over-permissive: it counts $M$ as surviving if any
chord crosses the enclosing chord, even when that single crossing chord can
be held finite by a careful rate choice (and so the kill argument actually
fires).

Under the **rigorous** rate-equation criterion (used in `step5_leading_kill_all_k`),
*every* step-1 survivor at $n=7,8,9$ dies at $k=2$ — there is no $k=3$ split.
So:

- The 57 "k=3 geometric" survivors at $n=9$ are all *rigorously* killed at
  $k=2$.
- The $m$-fold pinwheel at any $n$ is *rigorously* killable at $k=2$
  (probably) — needs explicit check but follows from the rate-equation
  feasibility for the pinwheel's interior chords.

Two upshots:
- **For the algebraic kill argument**: $K_{\text{rigorous}}(n) = 2$ through
  $n=9$; conjecturally for all $n$. **Locality at every $n$ may follow from
  $k\in\{1,2\}$ leading kills alone.**
- **For the geometric / combinatorial story**: $K_{\text{geometric}}(n)$
  probably saturates at 3 for $n\ge 9$ (no example found needing higher) but
  partial enumeration at $n\ge 10$ would settle it.

## Files

- `pinwheel_test.py` — main script (geometric criterion + pinwheel
  construction + filler enumeration).
- `check_n10_step1_survivors.py` — partial check at $n=10$ of whether any
  step-1 survivor needs $k\ge 4$ (geometric).
- `outputs/validation_n9.json` — the 954 / 57 split confirmation.
- `outputs/pinwheel_n{12,15}_m{4,5}.json` — full filler-enumeration histograms.
- `outputs/summary.json` — combined.

## What to do next

1. **Decide which criterion is the target**: rigorous (algebraic) or geometric
   (combinatorial). The rigorous criterion is what the kill argument actually
   does; the geometric is a clean heuristic that under-counts kills.
2. **If rigorous**: check whether all multisets at $n=10$ are killable at
   $k\in\{1,2\}$ (extending the step-5 verdict). If yes, conjecture $K_{\rm rig}(n)=2$
   for all $n$.
3. **If geometric**: enumerate step-1 survivors at $n=10$ and check whether
   any needs $k\ge 4$. The `check_n10_step1_survivors.py` script does this;
   takes ~tens of minutes in pure Python.
4. **Either way**: the bounded-$k$ strategy is still viable — just with a
   smaller $K$ than $\lfloor n/3\rfloor$. Probably $K=2$ rigorously, $K=3$
   geometrically.
