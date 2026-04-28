# n=9 Step-2 flip-graph connectivity check

**Backs paper §3 (Step-2 unitarity engine) and Item 7 of the master plan.**

This folder contains the verifier for whether Step-2's bare/special swap
identities, applied to the 9-gon's triangulations, connect them into a
single cyclic-orbit class — the property that gives unitarity at $n=9$.

## Headline result

| | Value |
|---|---:|
| Triangulations of the 9-gon | 429 = $C_7$ ✓ |
| Cyclic orbits | **49** (47 of size 9, 2 of size 3) |
| Step-2 swap edges (triangulation level) | 378 |
| Step-2 swap edges (orbit level) | 41 |
| **Connected components (orbit level)** | **16** |
| Largest component | 4 orbits |

> **Step-2 alone does NOT bridge all triangulation orbits at $n=9$.**
> The flip graph has 16 connected components: 9 of size 4, 6 of size 2,
> 1 of size 1.

The user's prompt expected 32 orbits and 1 component. The 32 came from a
different count — the 32 triangulation **external columns** that appear
in the n=9 cluster-partial matrix, which is a SUBSET of all 49 orbits.
The actual triangulation-orbit count under Z_9 is 49.

## Interpretation

### Step-2 is not sufficient at $n=9$

The bare/special swap is a single-edge move that exchanges
$X_{r, r+2} \leftrightarrow X_{r+1, r-1}$ (which always cross, so a
triangulation contains at most one of the two). Two triangulations are
swap-related iff they differ by exactly one such pair at some zone. At
$n \le 8$ the resulting flip graph is connected. **At $n=9$ it has
16 components.**

The 4 size-2-or-1 components (smaller than the typical 4-orbit size)
are likely the rotation-stabilised triangulations: with orbit sizes 3
and 9, certain swap moves that would normally produce length-4 cycles
collapse onto themselves.

### How the n=9 cluster matrix relates to this

The 90-orbit cluster analysis at
`../../step4_laurent_block_analysis/n9/` produced a 506-row
fingerprint matrix with 32 triangulation external columns (out of the
49 orbits). Cross-checking which of the 16 Step-2 components those 32
external triangulation orbits land in:

| | Count |
|---|---:|
| Step-2 components touched by ≥1 cluster external column | **12** |
| Components NOT touched (cluster has no relation involving them) | **4** |

So the cluster matrix can — in principle — provide *additional*
bridging equations on the 12 touched components, beyond what Step-2
alone gives. But **4 components are completely outside the cluster
matrix's reach**, so they need a separate mechanism (likely some
deeper Laurent / cluster equation, or a connection to the 21
non-cluster step-1 survivors).

### What this means for the n=9 locality theorem

The story is now:

| n=9 component | Mechanism | Status |
|---|---|---|
| Cluster orbits (90 of 113) | depth-1 fingerprint matrix, rank 90 = full | ✓ done |
| Non-cluster non-tri survivors (23) | single-orbit depth-1 cascade | ✓ done |
| **Triangulation orbits (49)** | **Step-2 (16 components) + cluster bridging on 12 of them** | **partially open** |

The "missing" piece is bridging the remaining components into a single
unitarity class. Possible avenues:
- The cluster matrix rows that involve the 4 untouched components
  may exist via deeper-Laurent fingerprints (depth-2 or higher).
- Higher-Laurent Step-2 identities (multi-bare swaps) may merge
  components.
- An "anchor" orbit at higher depth that simultaneously equates two
  different components.

This is the central open question for the n=9 closure.

## Files

- **`flip_graph_n9.py`** — verifier script (heavily commented per repo standards).
- **`outputs/n9_triangulation_orbits.json`** — manifest of the 49 cyclic orbits with reps.
- **`outputs/n9_step2_bare_swap_pairs.json`** — every Step-2 swap edge at the
  triangulation level and the orbit level, plus self-loops.
- **`outputs/n9_flip_graph_connectivity.txt`** — human-readable verdict, including
  per-component representatives.

## Run

```bash
python3 flip_graph_n9.py
```

Runs in ~1 second on a laptop.

## Cross-references

- Cluster analysis at $n=9$:
  `../../step4_laurent_block_analysis/n9/cluster_findings.md`
- Step-2 conservative classes at $n=7$ (sister folder):
  `../equivalence_classes/`
- Master plan Item 7 (unitarity from Step-2 flip-graph):
  root README §6.
