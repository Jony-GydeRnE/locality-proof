# `paper/` — formal write-ups

The core paper does **not yet exist** as a single document. The proof is
currently distributed across self-contained satellite write-ups in this
folder, each one targeting one section of the upcoming core paper. Once
the central conjecture (Item 2 in the root README) is
closed, these satellites will be folded into a single document.

For the precise theorem and proof strategy, read the root README first.

## The satellites, in paper-section order

| Subfolder | Backs paper § | What it is |
|---|---|---|
| [`elementary-geometric-background-and-understanding/`](elementary-geometric-background-and-understanding/) | §1. Introduction | Geometric / closed-timelike-curve story for newcomers — polygon ↔ Feynman-diagram dictionary, what a hidden zero is, why locality has to hold. |
| [`B=0->B_i = 0/`](B=0->B_i%20=%200/) | §3. Foundational lemma | If $B\|_{\mathcal Z_r}=0$, each homogeneity weight-component $B_i$ vanishes too. Generalizes Rodina's eqs ~15–17. Needed to apply the kill move term-by-term. |
| [`ex- proving locality&unitarity by hand/`](ex-%20proving%20locality%26unitarity%20by%20hand/) | §4. Worked examples | The current "core example" at $n=4, 5, 6$. Step 1 and Step 2 in full detail. **Doubles as the standing placeholder for the core paper** until it is written. |
| [`Laurent series for hard kills/`](Laurent%20series%20for%20hard%20kills/) | §6. Step-3 mechanism | The Laurent cascade kill at $n=7$. Proves all seven non-triangulation "fish" that escape Steps 1 and 2 die at the next Laurent order via a depth-1 cascade. |

## How the satellites map to the core paper

The intended core-paper outline is:

1. **§1. Introduction** — state the theorem (root README §1).
2. **§2. Setup** — polygon, chords, ansatz, mesh identity, hidden 1-zeros.
3. **§3. Foundational lemma** — `B=0->B_i = 0/`.
4. **§4. The kill mechanism** — Steps 1, 2, 3 statement + worked examples
   from `ex- proving locality&unitarity by hand/`.
5. **§5. Numerical evidence** — see `../computations/`.
6. **§6. The Laurent cascade** — `Laurent series for hard kills/`.
7. **§7. The combinatorial lemma** (Item 2, central open) — proof or
   strongest reduction.
8. **§8. The conditional locality theorem** — given Item 2, survivors
   $= C_{n-2}$ and the kill mechanism collapses the surviving space to
   one dimension spanned by $A_n^{\text{tree}}$.
9. **§9. Asymptotic locality** — survivor fraction $\to 0$ fast as $n\to\infty$.
10. **§10. Discussion** — extensions to massive case, gauge theories, loops.

Each subfolder is short enough to be a chapter, but self-contained
enough to read alone.

## Naming conventions

- Subfolder names are descriptive, not section-numbered, so adding new
  satellites later does not force renames.
- Each subfolder has its own README pointing into the file(s) and to the
  paper § it backs.
- Older drafts and superseded approaches live in
  `../old stuff/paper-older-works/`, never in here.
