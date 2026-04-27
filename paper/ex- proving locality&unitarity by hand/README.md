# §4 — Worked examples (and current core-paper placeholder)

**Backs paper §4 (worked examples) and Item 8 of the master plan.**

This is the most explicit by-hand walk-through of the kill mechanism in
the repo. It works through Steps 1 and 2 in full detail at
$n = 4, 5, 6$, with every coefficient identity written out. No Laurent
expansion is needed at these orders — that becomes essential only at
$n \ge 7$ (see `../Laurent series for hard kills/`).

**Until the core paper is written, this file is the standing placeholder
for it.** It is the closest thing to a self-contained proof at small $n$
that exists in the repo today.

## Files

- **`locality_unitarity_v5.tex`** — the worked example (April 2026, v5).
- **`locality_unitarity_v5.pdf`** — compiled (the file to read).

## What it covers

- §1. Setup. Polygon, chords, ansatz, mesh identity, hidden 1-zeros,
  master substitution, special chord $X_*$, free set $F_r$.
- §2. The kill move (statement + proof of Lemma).
- §3. $n = 4$ in full.
- §4. $n = 5$ in full: 15-dim ansatz, all 5 zones, 10 non-locals killed
  by Step 1 at three zones each, 5 unitarity equalities forming a
  5-cycle on the fan triangulations.
- §5. $n = 6$ in detail: 165-dim ansatz, 129/151 non-locals killed by
  Step 1, the remaining 22 falling into 5 cyclic orbits handled by
  cyclic-shifted Step 2, all 14 triangulation coefficients collapsed
  via flip-graph connectivity.

The file ends with a "general-$n$ outlook" stating the open conjectures
that the upcoming core paper will address.
