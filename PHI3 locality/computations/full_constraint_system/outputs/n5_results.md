# Full constraint-system rank at n = 5

- T(n) (multisets of size 2): **15**
- C_{n-2} (triangulations): **5**
- (k, r) pairs sampled: **10** (k = 1..2, r = 1..5)
- Constraints sampled: **14**
- **Rank** of constraint matrix: **14**
- **Nullspace dimension**: **1**
- Triangulation-indicator vector in nullspace: **True**
- Prime: 998244353
- Elapsed: 0.0 s

## Verdict

**Cokernel is 1-dim, spanned by the triangulation indicator.**

Every non-triangulation coefficient a_M is forced to zero, and
every triangulation coefficient equals one common scalar.
This is locality + unitarity at this n from the FULL k-zero
system -- strictly stronger than what Part I proves (k=1 only).
