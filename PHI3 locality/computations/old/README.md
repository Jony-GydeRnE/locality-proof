# `old/` — superseded experiments

These scripts supported the **older 15-dimensional ansatz / BCFW-induction
approach** to the proof, which has since been replaced by the
multiset-ansatz + Laurent direction (Steps 0–3 in the parent
`computations/README.md`).

They are kept here for provenance and because some observations they
produced may be cited as historical remarks. Do not treat them as part
of the active proof toolchain.

## Folders

| Folder | What it did | Why it is here, not active |
|---|---|---|
| [`independence_check/`](independence_check/) | Counts how many of the $n$ cyclic 1-zero loci are independent as constraints on the rational ansatz. Found $n-1$ at $n=6$ (correcting an earlier $\lfloor n/2 \rfloor$ conjecture). | Tied to the older 15-dim formulation; the new direction does not depend on a minimal independent set. |
| [`verify_factorization/`](verify_factorization/) | Numerical check that the BCFW residues at $n=5, 6$ factorize correctly. Supported the Part-II induction step in the old `main.tex`. | The current direction does not use BCFW induction at all. |
