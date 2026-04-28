# `paper/` — formal write-ups

The core paper does **not yet exist** as a single document. The proof is
currently distributed across self-contained satellite write-ups in this
folder. Each one targets one section / step of the upcoming core paper.

For the precise theorem and proof strategy, read the root README first.

## Top-level layout

```
paper/
├── README.md                                         ← you are here
├── elementary-geometric-background-and-understanding/  intro for newcomers
├── step1 Kill technique and statistics/               §4 + §5  step-1 layer-0 kill
├── step2 Laurent series for hard kills/               §6  Laurent cascade (n=7)
├── step3 block rule/                                  §7  block-rule kill (n=9)
└── Details missing in Hidden zero -> unitarity Rodina proof/
                                                       gap-filling lemmas in Rodina's setup
```

## The satellites, in paper-section order

| Subfolder | Backs paper § | What it is |
|---|---|---|
| [`elementary-geometric-background-and-understanding/`](elementary-geometric-background-and-understanding/) | §1. Introduction | Geometric / closed-timelike-curve story for newcomers. Polygon ↔ Feynman-diagram dictionary, what a hidden zero is, why locality has to hold. |
| [`step1 Kill technique and statistics/`](step1%20Kill%20technique%20and%20statistics/) | §4 + §5. Step-1 layer-0 kill | The original Step-1 kill mechanism: worked examples at $n=4,5,6$ (`locality_unitarity_v5`), the local-survival lemma proving Step-1 never kills a triangulation (`step-one-doesnt-kill-triangulations/`), and the global kill-rate statistics through $n=9$ (`step-one-kill-statistics/`). |
| [`step2 Laurent series for hard kills/`](step2%20Laurent%20series%20for%20hard%20kills/) | §6. Step-2 Laurent cascade | The depth-1 Laurent cascade kill that closes the seven $n=7$ "fish": one Laurent step deeper at one zone, with prerequisites all of layer-0 type at one neighbouring zone. |
| [`step3 block rule/`](step3%20block%20rule/) | §7. Step-3 block-rule kill | The block-rule cascade kill at $n=9$: when singleton cascade fails, the depth-1 fingerprint matrix on each connected component of the depth-1 cousin graph has trivial cluster nullspace. Verified at $n=9$ on a 90-orbit cluster (rank 90 = full). |
| [`Details missing in Hidden zero -> unitarity Rodina proof/`](Details%20missing%20in%20Hidden%20zero%20-%3E%20unitarity%20Rodina%20proof/) | §3. Foundational lemmas | Self-contained proofs of two lemmas Rodina either left as exercises or used implicitly: $B=0 \Rightarrow B_i=0$ (homogeneity decomposition) and the d-subset uniqueness of the master substitution. |

## How the satellites map to the core paper

The intended core-paper outline is:

1. **§1. Introduction** — state the theorem (root README §1).
2. **§2. Setup** — polygon, chords, ansatz, mesh identity, hidden 1-zeros.
3. **§3. Foundational lemmas** —
   `Details missing.../Proof of Rodina claim B=0->B_i = 0/` and
   `Details missing.../details for Rodina D subset argument/`.
4. **§4. The kill mechanism** — Step-1 layer-0 kill from `step1/`.
5. **§5. Numerical evidence + statistics** — `step1/step-one-kill-statistics/`
   and `../computations/`.
6. **§6. The Laurent cascade (Step-2)** — `step2 Laurent series for hard kills/`.
7. **§7. The block-rule kill (Step-3)** — `step3 block rule/`.
8. **§8. Conditional locality theorem** — given the block-rule kill at every $n$,
   survivors collapse to one dimension spanned by $A_n^{\text{tree}}$.
9. **§9. Asymptotic locality** — survivor fraction $\to 0$ fast as $n\to\infty$.
10. **§10. Discussion** — extensions to massive case, gauge theories, loops.

Each subfolder is short enough to be a chapter, but self-contained
enough to read alone.

## Naming conventions

- Subfolder names start with `step1 / step2 / step3 / Details missing / elementary`
  matching the proof's logical progression.
- Each subfolder has its own README pointing into the file(s) and to the
  paper § it backs.
- Older drafts and superseded approaches live in
  `../old stuff/paper-older-works/`, never in here.
