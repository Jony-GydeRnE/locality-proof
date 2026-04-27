# Step 0 — Setup sanity check

**Backs paper §2.** Before we try to prove the *converse*
($B|_{\mathcal Z_r} = 0$ for all $r$ ⟹ $B = c\cdot A_n^{\text{tree}}$),
we sanity-check the *forward* direction: does the Tr(φ³) tree amplitude
itself actually vanish on each cyclic 1-zero locus?

This is the easy direction of Rodina's hidden-zero theorem. Verifying it
numerically is a cheap baseline that catches setup bugs before any of the
later steps run.

## Subfolders

| Folder | What it does |
|---|---|
| [`verify_tree_zeros/`](verify_tree_zeros/) | Symbolically constructs $A_n^{\text{tree}}$ as a sum over triangulations and checks $A_n^{\text{tree}}\big|_{\mathcal Z_r} \equiv 0$ for every cyclic $r$. Runs in under a second. |
