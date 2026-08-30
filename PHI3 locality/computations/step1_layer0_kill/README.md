# Step 1 — Layer-0 kill mechanism

**Backs paper §4 + §5.** At each cyclic zone $\mathcal Z_r$, sending the
special chord $X_* = X_{r,r+2}$ and every chord outside $F_r$ to $\infty$
forces every coefficient $a_M$ with $M \subseteq F_r$ to vanish individually.
This is the simplest and most powerful kill in the proof — at $n=9$ it
already accounts for **99.888%** of all non-triangulation multisets.

The two experiments here characterize the mechanism *quantitatively*: how
many multisets does it kill, and how many escape?

## Subfolders

| Folder | What it does |
|---|---|
| [`kill_enumeration/`](kill_enumeration/) | For each $n$, enumerates every size-$(n-3)$ chord multiset and checks at which cyclic zones (if any) it is killed by Step 1. Reports the survivors — the inputs to Step 2 and Step 3. Headline: 0 non-tri survivors at $n \le 6$, 7 at $n=7$, 100 at $n=8$, 1 011 at $n=9$. |
| [`dual_X13_never_special/`](dual_X13_never_special/) | Repeats the enumeration with one zone *removed* (default $\mathcal Z_{1,3}$, so $X_{13}$ is never the special chord). The "extras" — multisets uniquely killed by the excluded zone — quantify each zone's individual contribution and are all non-triangulations, supporting the cyclic-symmetry expectation. |
