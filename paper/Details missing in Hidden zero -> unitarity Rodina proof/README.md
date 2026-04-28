# Details missing in Rodina's hidden-zero → unitarity argument

This folder collects self-contained proofs of two foundational lemmas
that Rodina (arXiv:2406.04234) either left as exercises or used
implicitly in the argument from the cyclic 1-zero conditions to
locality + unitarity. Filling these in is part of the master plan
(Item 1) and lives in §3 of the upcoming core paper.

## What's inside

| Subfolder | Lemma |
|---|---|
| [`Proof of Rodina claim B=0->B_i = 0/`](Proof%20of%20Rodina%20claim%20B%3D0->B_i%20=%200/) | **Homogeneity decomposition.** If $B$ vanishes on the cyclic 1-zero locus $\mathcal Z_r$, then so does each homogeneous weight-component $B_i$ separately. Generalises Rodina's eqs ~15–17 from the tree amplitude $A_n^{\text{tree}}$ to a fully general rational ansatz $B$. |
| [`details for Rodina D subset argument/`](details%20for%20Rodina%20D%20subset%20argument/) | **D-subset uniqueness.** The master substitution at zone $\mathcal Z_{r, r+2}$ is uniquely determined by the d-dimensional subset structure used to derive it; the argument is sketched in Rodina's paper but not fully written. |

Both lemmas are *needed* to apply the kill move term-by-term in the
ansatz framework (otherwise we couldn't equate weight-components or
reduce coefficients independently).

## Backs paper §3

The core paper §3 ("Foundational lemmas") imports both proofs.
