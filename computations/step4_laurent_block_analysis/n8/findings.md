# n=8 per-orbit audit — headline findings

This file summarises the key takeaways from `audit_orbits.py` /
`orbit_audit.json` / `orbit_audit_summary.md`, which audited *every*
zone × substitute pair $(Z, U)$ for each of the 13 cyclic-orbit
representatives, recording where a depth-1 Laurent cascade succeeds
(i.e., produces an equation whose every cousin is step-1 killable).

## 1. Depth-1 cascade is universal at $n=8$

> **Every one of the 13 orbits has at least one valid $(Z, U)$ recipe.**

| Orbit | size | # successful $(Z, U)$ depth-1 recipes |
|---:|---:|---:|
| 1 | 8 | 2 |
| 2 | 8 | 2 |
| 3 | 8 | 2 |
| 4 | 8 | 1 |
| 5 | 8 | 2 |
| 6 | 8 | 4 |
| 7 | 8 | 1 |
| 8 | 8 | 1 |
| 9 | 8 | 2 |
| 10 | 8 | 2 |
| 11 | 8 | 1 |
| 12 | 4 | 4 |
| 13 | 8 | 2 |
| **Total** | — | **26** |

Many orbits have multiple valid recipes (orbits 6 and 12 have four
each), giving the would-be all-$n$ proof slack: it suffices to
choose one zone per orbit, but the existence of multiple options
suggests the rule is robust.

## 2. Φ-I (the frame-configuration rule) is too restrictive

The Φ-I rule predicts a kill zone $\mathcal Z_{r,r+2}$ from the
missed-vertex set: it requires both endpoints $\{r, r{+}2\}$ to lie
in $V_{\text{missed}}(M)$. In the audit, **orbits 4, 7, 11, and 12 have
NO Φ-I frame candidates** (their missed-vertex set has no
cyclic-distance-2 pair), yet **each of these orbits still has at least
one valid depth-1 recipe**:

- Orbit 4 (V_missed = {2, 5}): recipe at $\mathcal Z_{5,7}$ — special $(5,7)$ has only vertex 5 in V_missed.
- Orbit 7 (V_missed = {5, 8}): recipe at $\mathcal Z_{5,7}$ — same situation.
- Orbit 11 (V_missed = {8}): recipe at $\mathcal Z_{8,2}$ — only one missed vertex.
- Orbit 12 (V_missed = {4, 8}): four recipes, all at zones whose special touches V(M).

**Conclusion.** A successful depth-1 cascade does *not* require the
zone's special chord to lie in $V_{\text{missed}}(M)$. The actually
load-bearing structural conditions are simpler:

> **Audit-derived rule.** A $(Z, U)$ pair gives a depth-1 cascade
> recipe if and only if:
>
> 1. $U \in M$ is a *non-bare substitute* on $Z$ (i.e., $U = X_{r+1, k}$ for some $k$),
> 2. its companion $Y_U = X_{r, k}$ is **not** in $M$ (K2 is clean for this pair),
> 3. the resulting depth-$1$ fingerprint equation's cousins are all step-1 killable.
>
> Condition (3) is empirical: across all $13 \times 8$ zone scans
> performed by the audit, *every* $(Z, U)$ satisfying (1) and (2)
> with $\ell_Z \ge 2$ produced a fingerprint equation; the question
> is just whether the cousins die — and at $n=8$, **they always do
> at depth 1** (sometimes via Φ-I, often not).

This is a strict generalisation of Φ-I; it should be the rule we
state in the paper as the conjectured all-$n$ recipe (Item 10).

## 3. Single Z_2-stabilised orbit (orbit 12) and the rest free

12 of the 13 orbits have size 8 (the action of $\mathbb Z_8$ is free).
Orbit 12, with representative $\{(1,3),(1,7),(2,6),(3,5),(5,7)\}$, has
size 4 — the $v \mapsto v + 4 \pmod 8$ shift fixes it. This orbit also
has the most recipes (four), one for every "pair-symmetric" zone
$\mathcal Z_{2,4}, \mathcal Z_{4,6}, \mathcal Z_{6,8}, \mathcal Z_{8,2}$.

## 4. ell_Z constraints

Across the 13 orbits, all valid recipes are at zones with $\ell_Z \ge 2$
(at least one bare or substitute, plus at least one more bare /
substitute / special). Zones with $\ell_Z \le 1$ never yield a recipe
(the Laurent expansion at sub-leading is degenerate).

**Specifically, the successful zones in the audit always have:**
- $n_{\text{bare}} \ge 1$ AND $n_{\text{nonbare\_subs}} \ge 1$, OR
- $n_{\text{nonbare\_subs}} \ge 2$.

i.e., the leading Laurent order is $\ge 2$, ensuring a non-trivial
sub-leading order with companion variables.

## 5. Where the data goes next

- These observations refine Item 10 (per-orbit depth-1 cascade
  uniformity, root README) into a sharper conjecture: the audit-derived
  rule above.
- The same audit at $n=9$ (113 orbits) is the next experiment;
  see `../cascade_n9/`.
- The full per-orbit JSON (`orbit_audit.json`) and per-survivor data
  (`recipe_analysis.md`) provide the input for any combinatorial proof
  of the rule.
