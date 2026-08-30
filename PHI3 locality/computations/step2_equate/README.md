# Step 2 — Bare↔special equate (the unitarity engine)

**Backs paper §4 + §5.** Step 1 *kills* coefficients; Step 2 *equates*
them. At each zone $\mathcal Z_r$ the bare relation
$X_{r+1,r-1} = -X_*$ glues pairs of multisets that differ only by this
swap, producing equalities $a_M = a_{M'}$. Cyclic shifts of these
equalities form a flip-graph on the $C_{n-2}$ triangulations; flip
connectivity of the associahedron $K_{n-1}$ then collapses every
triangulation coefficient to one common value.

The experiment here computes the equivalence classes induced by the
swap, both over-permissive and conservative versions, so we can see
exactly which equalities Step 2 produces directly versus which require
indirect chains.

## Subfolders

| Folder | What it does |
|---|---|
| [`equivalence_classes/`](equivalence_classes/) | For each $n$, partitions Step-1 survivors into bare↔special equivalence classes in two modes: PERMISSIVE (every swap, may over-equate) and CONSERVATIVE (only swap edges where both endpoints are survivors). Reports cyclic orbit-shapes and class sizes. The conservative mode gives full unitarity at $n=5$ and the lower-bound structure at $n \ge 6$. |
