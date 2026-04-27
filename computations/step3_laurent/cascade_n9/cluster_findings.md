# n=9 cluster analysis — findings (block-rule kill verified)

## Headline result

> **The 90-orbit coupled-cascade cluster around orbit 22 has cluster
> rank = 90 (full) at depth 1, residual r = 0.  All cluster orbits
> are forced to zero by the depth-1 fingerprint system alone, once
> external (triangulation + non-cluster-survivor) columns are
> separated.**

This is the *block-rule kill* the user asked to verify: orbit 22 (and
the other 3 cascade-failures 46, 88, 108) cannot be killed in isolation
at depth-1, but they are killed *collectively* by the depth-1 fingerprint
system on the cluster.

## Setup recap

- Per-orbit depth-1 cascade verified 94 / 113 of the n=9 cyclic-orbit
  reps. 4 are genuine failures (22, 46, 88, 108) — they have no depth-1
  recipe whose every cousin is step-1-killable.
- 15 are timeouts, of which spot-checks (orbits 28, 34, 35) showed
  depth-1 *does* close them given more time; orbit 25 is a genuine
  failure structurally identical to orbit 22.
- BFS from orbit 22 (a step-1 survivor whose depth-1 cousins are
  themselves step-1 survivors in *other* orbits) discovered a coupled
  cluster of 90 orbits before being terminated for the rank analysis.

## Matrix dimensions

| Quantity | Value |
|---|---|
| Cluster orbits $m$ | 90 |
| Depth-1 equations | 506 |
| External columns | 53 |
| Total columns | 143 |

External columns break down as:

| Kind | Count |
|---|---:|
| Triangulation cousins (tree-amplitude unknowns) | **32** |
| Non-cluster non-tri step-1 survivors | **21** |

## Rank / nullity

| Sub-matrix | Rank | Nullity |
|---|---:|---:|
| Full $[A \mid B]$ (cluster + external columns) | 140 | 3 |
| Cluster-only $A$ | **90** | **0** |

**Cluster residual $r = m - \mathrm{rank}(A) = 0$.**

Setting all external columns to zero (i.e., assuming the 21 non-cluster
survivors die via their own cascades, and triangulations are handled
by Step 2's flip-graph equality) forces **every** cluster orbit's
coefficient to vanish. The depth-1 system on the cluster alone is
already determining.

## What the residual 3 in the full system means

The 3-dimensional nullspace of $[A \mid B]$ corresponds to:

1. **Tree-amplitude direction** — the 32 triangulation cousins together
   carry the $A_n^{\text{tree}}$ coefficient, which is a single
   uncategorised scalar (Step 2's flip-graph kicks in here, not Step 1).
2. **Two non-cluster-survivor directions** — the 21 outside-cluster
   step-1 survivors that BFS would have eventually pulled in if it had
   continued. Once those are added, expect cluster nullity to still
   stay at 0 (or pick up at most the tree direction, depending on how
   triangulation columns combine).

The **$r = 0$ residual** is the key takeaway: the cluster as a whole is
fully determined by the depth-1 system, so the 4 failure orbits are not
genuinely uncaught — they're tied to the rest of the cluster through
linear constraints that, when solved jointly, force every $a_M = 0$.

## Comparison to n=6 perfect-matching cluster

At n=6 the four perfect matchings form a coupled cluster with
$\mathrm{rank}(M_{FP}) = 3$ and a single anchor (one Laurent step deeper
at one zone) closes the residual direction. At n=9, the analogous
cluster around orbit 22 is **dramatically larger** (90 orbits vs 4)
but **structurally simpler in one sense**: $r = 0$ already at depth-1
with no anchor needed (because the cluster's many depth-1 equations
overdetermine each other).

## Implications for Item 10

The conjecture (root README, Item 10) was:
> *every Step-1 survivor dies via a depth-1 Laurent cascade with all
> prerequisites of layer-0 type at one neighbouring zone.*

The n=9 data refines this:

> *the depth-1 cascade does not always isolate a single orbit; instead,
> for some orbits (like orbit 22 at n=9), the depth-1 fingerprint
> system collectively forces a coupled CLUSTER of orbits to zero.
> The block-rule kill at depth-1 closes the cluster jointly, even when
> no individual orbit-by-orbit recipe exists.*

This is consistent with the n=6 perfect-matching observation, but
generalises naturally to large clusters. The all-$n$ proof of locality
should now be framed as:
1. Step-1 kills coefficients with $M \in \mathcal K_{r,r+2}$ at any zone.
2. Step-2 equates triangulation coefficients via flip-graph connectivity.
3. **Cluster step-3**: groups of step-1 survivors form coupled clusters;
   the depth-1 fingerprint matrix on each cluster has trivial nullspace
   (modulo external columns which are themselves killed by other clusters
   or by step-2).
4. Triangulation coefficients collapse to one common scalar (= the tree).

## Next steps

- Continue the BFS to absorb the 21 non-cluster survivors and verify
  the cluster keeps $r = 0$.
- Check whether the cluster covers ALL 113 step-1-survivor orbits
  (likely, given the growth pattern).
- Once full BFS converges, compare to the full $n=9$ symbolic nullspace
  (= the equivalent of `full_ansatz_nullspace/full_ansatz_n7.py` but for
  $n=9$) — they should agree up to the triangulation tree direction.
