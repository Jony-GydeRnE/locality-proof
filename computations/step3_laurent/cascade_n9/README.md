# Cyclic-orbit decomposition + cascade-on-reps at $n=9$

**Companion to `../cascade_n8/`. Provides empirical evidence for Item 10
(per-orbit depth-1 cascade uniformity) at $n=9$.**

## Strategy

At $n=9$ there are 1 011 non-triangulation step-1 survivors. Running the
full Laurent cascade on all 1 011 would take ~3 hours of sympy work.
But the kill mechanism is cyclically equivariant: a depth-1 recipe
$(Z, k, U, Y_U, \mathrm{fp})$ for one orbit member $M$ shifts to a
recipe for any cyclic image $M + s$ at $(Z + s, k, U + s, Y_U + s,
\mathrm{fp} + s)$. So we only need one cascade run per
$\mathbb Z_9$-orbit.

This folder runs that strategy in two stages.

## Stage 1 — orbit decomposition (`orbit_decomposition.py`)

Enumerates the 1 011 step-1 survivors, partitions them into cyclic
orbits, and dumps one canonical representative per orbit.

> **Result: 113 cyclic orbits at $n=9$** — 112 of size 9 (the action
> is free) and 1 of size 3 (stabilised by $\langle v \mapsto v + 3 \rangle$).

| Orbit size | # orbits | survivor coverage |
|---|---:|---:|
| 9 | 112 | 1 008 |
| 3 | 1 | 3 |
| **Total** | **113** | **1 011** ✓ |

Outputs: `orbits_n9.json` (manifest for the cascade-on-reps stage)
and `orbits_n9.md` (human-readable summary).

## Stage 2 — cascade on orbit reps (`cascade_kill_orbit_reps_n9.py`)

For each of the 113 representatives, runs `search_depth1_recipe`
(imported from `../cascade_n8/cascade_kill_n8.py`, which is parametric
in $n$). Logs:

- orbit size, $M_{\text{rep}}$, $V_{\text{missed}}(M_{\text{rep}})$;
- recipe: kill zone $Z$, Laurent order $k$, substitute $U$, companion
  $Y_U$, fingerprint $X_{Y_U} / (\dots)$;
- whether $Y_U \in M_{\text{rep}}$ (K2-violation flag) at the chosen
  zone — almost always `False` if a depth-1 recipe is found.

By cyclic equivariance, success on the rep lifts to success on every
member of its orbit, so 113 successful cascades certify all 1 011
survivors.

Outputs: `results_cascade_n9_reps.txt` (human-readable per-rep trace)
and `results_cascade_n9_reps.json` (machine-readable record per rep).

### Status

- Smoke test (2 reps): **2 / 2 succeeded** in 5.6 min, average ~167 s/rep.
- Full run (113 reps): in progress; estimated ~5 hours total.

## Files

| File | Purpose |
|---|---|
| `orbit_decomposition.py` | Stage 1 — enumerate survivors, decompose into 113 orbits. |
| `orbits_n9.json` | Orbit manifest (one rep per orbit). |
| `orbits_n9.md` | Human-readable orbit summary. |
| `cascade_kill_orbit_reps_n9.py` | Stage 2 — run cascade on the 113 reps. |
| `results_cascade_n9_reps.txt` | Per-rep cascade trace (written incrementally). |
| `results_cascade_n9_reps.json` | Per-rep record with recipe details. |

## Run

```bash
python3 orbit_decomposition.py            # ~10 s; writes orbits_n9.json
python3 cascade_kill_orbit_reps_n9.py     # ~5 h; writes results_*
python3 cascade_kill_orbit_reps_n9.py 5   # smoke test on first 5 reps
```

## What this means for Item 10

If all 113 representatives admit a depth-1 cascade with K2-clean
substitutes (Y_U ∉ M_rep), the empirical case for the all-$n$
conjecture extends from $n=7$ (7/7) and $n=8$ (100/100) to $n=9$
(1011/1011), and the recipes give us the orbit-by-orbit data needed
to formulate a uniform combinatorial rule.
