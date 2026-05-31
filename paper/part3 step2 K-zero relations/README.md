# Part III — classification of Step-2 K-zero relations

`step2_part3.tex` → `step2_part3.pdf` (6 pp).

**Core point.** Vanishing of $B_n$ on a k-zero $\mathcal Z^{(k)}$ carries more information than the
Part-II asymptotic kill. Two graded decompositions apply:

1. **Boundary-propagator grading** (BCFW homogeneity, note 12): $B_n|_{\mathcal Z^{(k)}}=\sum_m B_n^{(m)}$
   with each $B_n^{(m)}$ vanishing independently.
2. **D-subset decomposition** (note 18): within each $B_n^{(m)}$, factor out the born-free monomial;
   each D-subset vanishes independently.

Then the **X→0 cuts** (the "Rodina cuts") read off explicit linear relations among the surviving
coefficients $a_M$. We classify these for $k=1,2,3$ at $m=1,2,3$.

**Universal relation** (Theorem 1, the binomial trace): for every $k,n,m\ge1$,
$$\sum_{j=0}^m (-1)^j a_{s^{m-j} b^j} = 0,$$
where $s=X_{1,k+2}$ (special), $b=X_{k+1,n}$ (bare). At $m=1$ this is the unitarity relation
$a_{X_{1,k+2}}=a_{X_{k+1,n}}$ (recovering $a_{13}=a_{2n}$ at $k=1$). At $m=2,3$ for $k\ge2$ these
are **genuinely new** relations, e.g.
$$a_{14,14}+a_{3n,3n}-a_{14,3n}=0,\quad
  a_{14,14,14}-a_{14,14,3n}+a_{14,3n,3n}-a_{3n,3n,3n}=0.$$

Plus per-column and per-offset relations (Theorem 2: structure at order m). The 2-zero, 3-zero are
worked out explicitly.

**Build:** `python3 figures/make_figures.py` then `tectonic step2_part3.tex`.

**Source notes:** 12 (BCFW grading), 18 (D-subsets), 26 (Step-2 for 1-zero), 27 (Step-2 for 2- and
3-zero — the load-bearing note for this paper).

Companion to Part II (`paper/part2 fat-zero generalization/`) — together they form "Paper 2" for
Rodina.
