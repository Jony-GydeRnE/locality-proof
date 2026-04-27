# Laurent cascade kill at $n = 8$ — depth-1 verification

**Companion to `paper/step-one-kill-statistics/` and `paper/Laurent series for hard kills/`.**
**Provides empirical support for Item 10 (the central open conjecture).**

## Headline result

> **All 100 of the $n=8$ Step-1 survivors die via a depth-1 Laurent cascade.**
> Zero failures. Verified symbolically with sympy in 4 536 s ($\approx 76$ min).

This is exactly the same phenomenon as the seven $n=7$ "fish" (proven
in `../cascade_n7/cascade_kill_n7.py`), now confirmed to persist at
$n=8$: every Step-1 survivor admits a kill recipe of the form

> *one Laurent step deeper at one zone, with all prerequisites of
> layer-0 type at one neighbouring zone.*

If this depth-1 phenomenon continues to hold at every $n \ge 7$
(Item 10 of the master plan), the all-$n$ locality theorem reduces to
verifying it on each cyclic orbit, which is a finite per-orbit
combinatorial check.

## What the script does (quick description)

For each of the 100 step-1 survivors at $n=8$, the script searches
over candidates `(zone Z, Laurent order k, fingerprint fp)` for a
triple satisfying:

1. The survivor $M$ contributes to the order-$k$ Laurent expansion of
   $1/\prod_{c\in M} X_c$ at zone $Z$ with nonzero scalar on the
   free-variable monomial `fp`.
2. Every *other* multiset that contributes to the same Laurent
   coefficient is itself step-1 killable at some (possibly different)
   zone.

If both hold, substituting the cousins' zeros into the Laurent
equation collapses it to `scalar(M) * a_M = 0`, hence $a_M = 0$.

The candidate generation mimics the worked $n=7$ recipe: pick a
substitute chord in $M$, take its companion as the fingerprint
numerator, and the rest of $M$'s free chords as the denominator. Try
$k = \text{leading order} + 1$ and $\text{leading order} + 2$.

## Distribution of kill recipes across the 100 survivors

| Laurent order $k$ | # survivors |
|---|---:|
| 3 (i.e., leading + 1) | 90 |
| 4 (i.e., leading + 2) | 10 |

By kill zone (rough cyclic symmetry):

| Zone | Survivors killed there |
|---|---:|
| $\mathcal Z_{1,3}$ | 25 |
| $\mathcal Z_{2,4}$ | 21 |
| $\mathcal Z_{3,5}$, $\mathcal Z_{4,6}$, $\mathcal Z_{5,7}$, $\mathcal Z_{6,8}$ | 10–11 each |
| $\mathcal Z_{1,7}$, $\mathcal Z_{2,8}$ | 5 each |

The asymmetry across zones is real and reflects the script's search
order (smallest $r$ first); other zones would give equally valid
recipes for the same survivors by cyclic rotation.

## Files

- **`cascade_kill_n8.py`** — the verifier (heavily commented; every
  function has a `LOGIC` and a `PHYSICS / MATHEMATICS` docstring per
  the repo contribution standards in the root README).
- **`results_cascade_n8.txt`** — full human-readable trace (100 blocks,
  one per survivor: kill zone, fingerprint, equation, cousin kills).

## Run

```bash
python3 cascade_kill_n8.py            # full run (~76 min)
python3 cascade_kill_n8.py 5          # smoke test on first 5 survivors
```

## What this means for the proof

- **Items 1, 2, 6, 8** of the master plan are closed.
- **Item 5** (numerical evidence at $n=7, 8$): closed for $n=7$ and
  $n=8$ — both completely verified.
- **Item 10** (per-orbit depth-1 cascade uniformity): **central open**;
  now empirically verified at $n=7$ (all 7 fish, hand-picked recipe)
  and $n=8$ (all 100 survivors, automated recipe search). $n=9$
  (1 011 survivors) is the next experiment.

If $n=9$ also shows uniform depth-1, the empirical case for Item 10
becomes very strong, and the remaining work is the combinatorial
proof that the depth-1 recipe is *always* available.
