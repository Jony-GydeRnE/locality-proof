# step 4 — Laurent block analysis (n=8 and n=9)

**Backs paper §6 + §7 (Laurent cascade and block-rule kill).
Provides empirical evidence for the all-$n$ locality conjecture (Item 10
of the master plan).**

n=7 was the original "single-orbit cascade" demonstration; that lives in
[`../step3_laurent/cascade_n7/`](../step3_laurent/cascade_n7/).
This folder takes the cascade analysis to n=8 and n=9, where the picture
gets richer.

## At a glance

| $n$ | Survivors | Method | Result |
|---|---:|---|---|
| 7 | 7 "fish" | depth-1 singleton cascade, hand-picked recipe | 7/7 die individually (`../step3_laurent/cascade_n7/`) |
| **8** | 100 | depth-1 singleton cascade, automated search | **100/100 die individually** ([`n8/`](n8/)) |
| **9** | 1 011 / 113 orbits | depth-1 singleton cascade per orbit rep | 94/113 die individually; 4 genuine failures, 15 timeouts |
| **9 cluster** | 90-orbit cluster around orbit 22 | **block-rule kill** via fingerprint matrix | rank 90 = m, **nullity 0** ([`n9/`](n9/)) |

The headline:

> **At $n=9$ the singleton-cascade rule *sometimes* fails at the
> orbit level, but the BLOCK-RULE Laurent kill — many depth-1
> fingerprint equations taken jointly — still forces every cluster
> orbit to zero.**

This is the central conceptual shift the work documents.

## Why the n=8 → n=9 transition matters

At $n=7, 8$ every step-1 survivor admits its own depth-1 cascade recipe:
one Laurent step deeper at one zone, with all prerequisites of layer-0
type at one neighbouring zone. The all-$n$ conjecture (Item 10) was
originally framed as "this singleton recipe persists to all $n$".

At $n=9$ that's no longer true orbit-by-orbit. Four orbits
(22, 46, 88, 108) have *no* depth-1 recipe whose every cousin is
step-1 killable — because the cousins are themselves step-1 survivors.
This is the "coupled subsystem" scenario.

But running BFS in the depth-1 cousin graph from orbit 22 collects a
**90-orbit cluster** (out of 113), and the depth-1 fingerprint matrix
on that cluster has **rank 90 = full**. With external columns set to 0
(triangulations and non-cluster survivors handled by their own
mechanisms), every cluster orbit is forced to zero. The four "failures"
are killed jointly, not individually.

The refined Item 10 reads:

> *The depth-1 fingerprint system, on each connected component of
> the depth-1 cousin graph, has trivial cluster nullspace.*

This is a strict generalisation of the singleton rule; singletons are
just clusters of size 1.

## Subfolders

| Folder | Contents |
|---|---|
| [`n8/`](n8/) | Singleton-cascade verification at n=8: 100/100 succeed at depth 1. Per-orbit audit (13 orbits, 26 valid recipes). Φ-I rule analysis (48/100 match). |
| [`n9/`](n9/) | Full cascade run on 113 orbit reps; deep diagnostic on anomalies; 90-orbit cluster matrix analysis with rank/nullity computation; block-rule kill verification. |

## Connection to n=6 perfect-matching cluster

The n=6 perfect-matching cluster has 4 orbits, $M_{FP}$ rank 3, and a
single anchor closes the residual direction. The n=9 cluster is
*much larger* (90 orbits) but *structurally simpler*: the depth-1
fingerprint system is so over-determined that $r = 0$ already, with
no anchor needed. The pedagogical pattern is the same — coupled
clusters whose collective nullspace dies — but the dimensionality
varies with $n$.
