# Step 3 — Laurent cascade for hard kills

**Backs paper §6.** Through $n = 6$, Steps 1 and 2 alone kill every
non-local coefficient. At $n = 7$, exactly seven non-triangulation
multisets — the cyclic family of "fish" — escape both. Step 3 catches
them by going to higher Laurent orders in $1/X_*$.

The mechanism: at zone $\mathcal Z_{r,r+2}$, expand
$\frac{1}{X_{r,k}-X_*} = -\frac{1}{X_*} - \frac{X_{r,k}}{X_*^2} - \frac{X_{r,k}^2}{X_*^3} - \cdots$
and read off coefficients of free-variable monomials at successive orders.
Each fish dies at order $X_*^{-3}$ (depth 1) via a 6-term fingerprint
equation whose other 5 contributors ("cousins") are all killed by Step 1
at one *neighbouring* zone. Substituting those zeros isolates the fish.

## Subfolders

| Folder | What it does |
|---|---|
| [`cascade_n7/`](cascade_n7/) | Symbolic verification of the depth-1 cascade kill for all seven $n=7$ fish. Reproduces the fingerprint equation, identifies the kill zone for each cousin, and confirms the equation collapses to $a_M = 0$ for every fish. Output trace in `results_cascade_n7.txt`. |

The same shape (one Laurent step deeper at one zone, prerequisites all
of Layer-0 type at one neighbouring zone) is conjectured to handle every
Step-1 survivor at $n \ge 8$ as well; this is the next experiment to
extend.
