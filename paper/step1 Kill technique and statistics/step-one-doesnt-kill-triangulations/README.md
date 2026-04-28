# §3.5 / §4 — Step 1 doesn't kill any triangulation (Item 2 closed)

**Backs paper §4 (kill mechanism — local survival of triangulations).
Closes Item 2 of the master plan.**

This satellite contains the proof that **no triangulation $T$ of the
$n$-gon is killable by the layer-0 (Step-1) hidden-zero mechanism at any
cyclic zone $\mathcal Z_{r,r+2}$.**

The argument is two pages, three exhaustive cases, and uses *only* the
fact that the polygon edge $(r, r{+}1)$ sits in exactly one triangle of
$T$. No Meisters two-ears theorem is needed.

## What this proves vs. what was originally conjectured

The original Item 2 in the master plan was the combinatorial
conjecture:

> Every $(n-3)$-element subset of the free set $F_r$ contains a
> crossing pair of chords  ⟺  no triangulation $T$ satisfies
> $T \subseteq F_r$.

The theorem proved here is **strictly stronger**:

> $T \notin \mathcal K_{r,r+2}$ for every triangulation $T$ and every
> $r$.

That is, $T$ fails (K1) or (K2) at every zone — ruling out not just
$T \subseteq F_r$, but also the case where $T$ contains a substitute
chord $X_{r+1,k}$ without its companion $X_{r,k}$. The original
conjecture is an immediate corollary: if $T \subseteq F_r$ then $T$
trivially satisfies (K1) and (K2), hence $T \in \mathcal K_{r,r+2}$,
contradicting the theorem.

## Files

- **`triangulation_layer0.tex`** — the proof.
- **`triangulation_layer0.pdf`** — compiled.

## The argument in one sentence

The polygon edge $(r, r{+}1)$ is a side of exactly one triangle of $T$;
let $a$ be its third vertex. Three cases on $a$ each force a different
violation:

| Case | $a$ | Violation |
|---|---|---|
| 1 | $a = r{+}2$ | side $(r, r{+}2)$ is the **special** chord and is in $T$ → (K1) fails |
| 2 | $a = r{-}1$ | side $(r{+}1, r{-}1)$ is the **bare** chord and is in $T$ → (K1) fails |
| 3 | $a \in \{r{+}3, \dots, r{+}n{-}2\}$ | sides $(r, a)$ and $(r{+}1, a)$ are the **(companion, substitute) pair** $(X_{r,a}, X_{r+1,a})$, both in $T$ → (K2) fails |

The three cases exhaust the $n-2$ possible values of $a$, so the
conclusion holds for every triangulation. ∎
