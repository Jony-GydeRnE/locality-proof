# Full constraint-system rank at n = 6

- T(n) (multisets of size 3): **165**
- C_{n-2} (triangulations): **14**
- (k, r) pairs sampled: **18** (k = 1..3, r = 1..6)
- Constraints sampled: **164**
- **Rank** of constraint matrix: **164**
- **Nullspace dimension**: **1**
- Triangulation-indicator vector in nullspace: **True**
- Prime: 998244353
- Elapsed: 0.2 s

## Verdict

**Cokernel is 1-dim, spanned by the triangulation indicator.**

Every non-triangulation coefficient a_M is forced to zero, and
every triangulation coefficient equals one common scalar.
This is locality + unitarity at this n from the FULL k-zero
system -- strictly stronger than what Part I proves (k=1 only).
