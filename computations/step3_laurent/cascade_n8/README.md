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

## Recipe analysis (`analyze_recipes.py`)

The follow-up script `analyze_recipes.py` parses the cascade trace
and computes per-survivor structural data:

- **Cyclic orbit decomposition.** The 100 survivors break into **13
  $\mathbb Z_8$-orbits**: 12 orbits of size 8 (free) and 1 orbit of
  size 4 (stabilised by the $\mathbb Z_2$ subgroup). All orbit sizes
  divide 8 as required.
- **Missed-vertex set $V_{\text{missed}}(M) = \{1,\dots,8\} \setminus V(M)$.**
- **Frame candidates** (cyclic-distance-2 pairs inside $V_{\text{missed}}$)
  — these are the kill-zone candidates predicted by the **Φ-I rule**:
  *the kill zone $\mathcal Z_{r,r+2}$ has $(r, r{+}2) \subseteq V_{\text{missed}}(M)$.*
- **Φ-I prediction match against the empirical kill zone.**

### Φ-I match: 48 / 100 (per-orbit breakdown)

| Orbit | Canonical rep | Size | Φ-I match |
|---:|---|---:|---:|
| 1 | $\{(1,3),(1,3),(1,7),(3,5),(5,7)\}$ | 8 | **8 / 8** |
| 2 | $\{(1,3),(1,4),(1,7),(3,5),(5,7)\}$ | 8 | **8 / 8** |
| 3 | $\{(1,3),(1,4),(1,7),(4,6),(5,7)\}$ | 8 | 6 / 8 |
| 4 | $\{(1,3),(1,4),(1,7),(4,6),(6,8)\}$ | 8 | 0 / 8 |
| 5 | $\{(1,3),(1,4),(2,8),(4,6),(6,8)\}$ | 8 | 6 / 8 |
| 6 | $\{(1,3),(1,4),(3,8),(4,6),(6,8)\}$ | 8 | 2 / 8 |
| 7 | $\{(1,3),(1,6),(1,7),(2,4),(4,6)\}$ | 8 | 0 / 8 |
| 8 | $\{(1,3),(1,6),(1,7),(3,5),(4,6)\}$ | 8 | 0 / 8 |
| 9 | $\{(1,3),(1,6),(1,7),(3,5),(5,7)\}$ | 8 | 6 / 8 |
| 10 | $\{(1,3),(1,7),(2,4),(3,5),(5,7)\}$ | 8 | 6 / 8 |
| 11 | $\{(1,3),(1,7),(2,4),(4,6),(5,7)\}$ | 8 | 0 / 8 |
| 12 | $\{(1,3),(1,7),(2,6),(3,5),(5,7)\}$ | 4 | 0 / 4 |
| 13 | $\{(1,3),(1,7),(3,5),(5,8),(6,8)\}$ | 8 | 6 / 8 |

### What the 48% match rate means

**Φ-I as currently stated is not the full structural rule.** It works
perfectly on orbits 1 and 2 (where the survivor's missed-vertex set
contains the chosen kill-zone special), partially on orbits 3, 5, 9,
10, 13 (cyclic action with multiple frame candidates per survivor;
the empirical recipe picks one but not always the Φ-I-predicted one),
and **never** on orbits 4, 7, 8, 11, 12 — these survivors are killed
at zones whose special chord has at least one endpoint in $V(M)$.

This empirical data:

1. **Refines the all-$n$ conjecture (Item 10).** A uniform recipe rule
   needs to handle the 0% orbits, where the kill zone touches $V(M)$.
   Φ-I missed those because its hypothesis (frame entirely in
   $V_{\text{missed}}$) doesn't hold there.
2. **Suggests a richer structural rule.** Possibly a refinement: pick
   a substitute chord $c_{\text{sub}} = (r{+}1, k) \in M$ and use its
   companion as the fingerprint, with kill zone $\mathcal Z_{r, r+2}$ —
   regardless of whether the special is in $V_{\text{missed}}$.
3. **Justifies the per-orbit verification approach.** Each of the 13
   orbits has a single recipe pattern (under cyclic shift); the
   all-$n$ proof needs only to verify one recipe per orbit.

## Per-orbit structural audit (`audit_orbits.py`)

`audit_orbits.py` does a deeper structural audit than `analyze_recipes.py`.
For each of the 13 orbit representatives, it walks every zone $Z$ and
every substitute $U \in M_{\text{rep}}$ at $Z$ with companion $Y_U \notin M$
(K2-clean), and tries the depth-1 cascade at $(Z, U)$.

> **Headline.** Every one of the 13 orbits has at least one valid
> depth-1 $(Z, U)$ recipe — total 26 valid recipes across all orbits.
> See `orbit_audit_findings.md` for takeaways. The audit confirms the
> Φ-I frame rule is too restrictive (4 orbits with empty Φ-I but with
> recipes anyway) and suggests a strictly weaker rule: every K2-clean
> $(Z, U)$ with $\ell_Z \ge 2$ produces a recipe whose cousins all
> step-1-die.

## Files

| File | What it is |
|---|---|
| `cascade_kill_n8.py` | Cascade verifier (one depth-1 recipe per survivor). |
| `results_cascade_n8.txt` | Full cascade trace (100 survivor blocks). |
| `analyze_recipes.py` | Per-survivor analyser: orbits, missed vertices, Φ-I match. |
| `recipe_analysis.md` | Per-survivor markdown table. |
| `orbits.md` | Orbit decomposition (13 orbits). |
| `analysis_summary.txt` | Short summary of the per-survivor analysis. |
| `audit_orbits.py` | Per-orbit audit: every (Z, U) recipe per orbit rep. |
| `orbit_audit.json` | Machine-readable audit record (zones × substitutes × recipes). |
| `orbit_audit_summary.md` | Per-orbit audit table. |
| `orbit_audit_findings.md` | **Headline takeaways from the audit.** |

## Run

```bash
python3 cascade_kill_n8.py            # full cascade run (~76 min)
python3 cascade_kill_n8.py 5          # smoke test on first 5 survivors
python3 analyze_recipes.py            # parse + analyse (seconds)
python3 audit_orbits.py               # per-orbit audit (~10 min)
```
