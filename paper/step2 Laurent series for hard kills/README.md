# §6 — The Laurent cascade for hard kills (Step-2 in the new naming)

**Backs paper §6.** The cascade kill that handles step-1 survivors at
$n = 7$ (and at $n = 8$ as a singleton kill). At $n = 9$ singleton
cascade fails; the **block-rule** generalisation lives in
[`../step3 block rule/`](../step3%20block%20rule/).

Through $n = 6$, Steps 1 and 2 of the kill mechanism alone kill every
non-local coefficient. At $n = 7$, exactly **seven** non-triangulation
multisets — the cyclic family of "fish" — escape both. This document
proves all seven die at the next Laurent order, via a depth-1 cascade
of remarkably simple shape:

> *one Laurent step deeper at one zone, with prerequisites all of
> Layer-0 type at one neighbouring zone.*

Concretely: at zone $\mathcal Z_{5,7}$ for the canonical fish $M_1$,
the order-$X_*^{-3}$ Laurent expansion of $1/\prod X_c$ produces a
6-term fingerprint equation. Five of the six contributors are killed
by Step 1 at $\mathcal Z_{4,6}$ — the *one* neighbouring zone — and
the equation collapses to $a_{M_1} = 0$. The other six fish are cyclic
rotations and inherit the same recipe.

The companion script in
[`../../computations/step3_laurent/cascade_n7/`](../../computations/step3_laurent/cascade_n7/)
verifies every step symbolically with sympy.

## Files

- **`cascade_n7.tex`** — the proof.
- **`cascade_n7.pdf`** — compiled (the file to read).

## Outlook

The paper closes by conjecturing that **every** Step-1 survivor at
$n \ge 8$ dies via a depth-1 cascade of the same shape. Verifying this
at $n = 8, 9$ is the next experiment to run; if it holds, the path to
an all-$n$ locality proof reduces to a single-step combinatorial
induction.
