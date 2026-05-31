# Part I — the 1-zero asymptotic analysis

`one_zero_part1.pdf` (7 pp).

**Core point.** Locality of $\mathrm{Tr}(\phi^3)$ tree amplitudes from the cyclic hidden
**1-zero** alone, proven through $n=9$.

Contents:
1. **Setup.** $n$-gon chords, ansatz $B_n=\sum_M a_M/\prod X_c$, hidden-zero surface
   $\mathcal Z_{1,3}$ (special $X_{13}$, bare $X_{2n}$), the choice function $\sigma$ and
   the free chord set $F_\sigma$.
2. **The kill theorem (Theorem 1).** Asymptotic limit $\Omega_\sigma$ sends constrained
   chords to $\infty$; survivors form a Laurent polynomial which must vanish, so
   $a_M=0$ for every $M\subseteq F_\sigma$. Conditions (K1)/(K2).
3. **Triangulations escape (Theorem 2).** The unique triangle containing leg $X_{12}$
   pins the special, the bare, or a constrained pair into every triangulation $T$.
4. **Asymptotic completeness.** Union bound gives $S(n)/T(n)\le 8/n+O(n^{-2})$.
   Empirical kill % through $n=9$: 100%, 100%, 99.70%, 99.76%, 99.89%.
5. **Step-1 survivors structure.** $(S1)$ at least one short chord, $(S2)$ every
   vertex witnessed (Proposition 2). The $n=7$ "fish" (Theorem 3): Laurent
   fingerprint at one zone whose non-target contributors are Step-1 killed at a
   neighbouring zone.
6. **Beyond $n=7$.** $n=8$: 100 survivors, one cyclic orbit, depth-1 cascade.
   $n=9$: 1011 survivors, 113 orbits, 23 single-orbit cascades + 90-orbit **block
   rule** (a $506\times 90$ Laurent-coefficient matrix of full rank 90).

**Status of locality through $n=9$:** proven. The all-$n$ conjecture is left open
in this paper.

**Companion:** Part II (`paper/part2 fat-zero generalization/`) extends the kill
mechanism from the 1-zero to the $k$-zero; Part III (`paper/part3 step2 K-zero
relations/`) classifies the Step-2 coefficient equations across the $k$-zeros.
