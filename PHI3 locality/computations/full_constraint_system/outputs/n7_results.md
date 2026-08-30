# Full constraint-system rank at n = 7

- T(n) (multisets of size 4): **2380**
- C_{n-2} (triangulations): **42**
- (k, r) pairs sampled: **28** (k = 1..4, r = 1..7)
- Constraints sampled: **2379**
- **Rank** of constraint matrix: **2379**
- **Nullspace dimension**: **1**
- Triangulation-indicator vector in nullspace: **True**
- Prime: 998244353
- Elapsed: 70.4 s

## Verdict

**Cokernel is 1-dim, spanned by the triangulation indicator.**

Every non-triangulation coefficient a_M is forced to zero, and
every triangulation coefficient equals one common scalar.
This is locality + unitarity at this n from the FULL k-zero
system -- strictly stronger than what Part I proves (k=1 only).
