# §4 + §5 — Step-1 (layer-0) kill: technique, local-survival lemma, statistics

This folder collects everything about the **Step-1 layer-0 hidden-zero
kill mechanism** — the simplest of the kill steps and the one that
handles the overwhelming majority of non-local coefficients at every
$n$ measured.

## What's inside

| File / subfolder | What it is |
|---|---|
| `locality_unitarity_v5.tex` / `.pdf` | **The current "core example" placeholder.** A by-hand walk-through of Steps 1 and 2 at $n=4, 5, 6$, with every coefficient identity written out. The tightest self-contained proof at small $n$. |
| [`step-one-doesnt-kill-triangulations/`](step-one-doesnt-kill-triangulations/) | **Local-survival lemma.** No triangulation $T$ is ever Step-1 killable at any zone. Two-page proof, three exhaustive cases, using only the polygon edge $(r, r{+}1)$. (Closes Item 2 of the master plan.) |
| [`step-one-kill-statistics/`](step-one-kill-statistics/) | **Global kill-rate data.** Counts of step-1 survivors through $n=9$ (99.888% non-tri kill at $n=9$, monotonically increasing rate) and the asymptotic conjecture $S(n)/T(n) \to 0$. |

## Backs paper §4 + §5

- **§4. The kill mechanism** — Step-1 layer-0 kill is the *first* of the
  three steps; the worked example file (v5) explains the move at small
  $n$, and the local-survival lemma proves the kill never accidentally
  catches a triangulation coefficient.
- **§5. Numerical evidence** — the kill-statistics subfolder gives the
  global kill-rate table (n = 4..9), the cyclic-orbit decomposition of
  survivors, and the asymptotic conjecture.

## Where the supporting computations live

- Step-1 enumeration scripts:
  `../../computations/step1_layer0_kill/kill_enumeration/`
- "Dual" experiment ($X_{13}$ never special):
  `../../computations/step1_layer0_kill/dual_X13_never_special/`
- Hand-curated survivor diagrams ($n=8, 9$):
  `../../computations/survivor_gallery/`

## Read order

1. `locality_unitarity_v5.pdf` — see Steps 1 and 2 in action at $n \le 6$.
2. `step-one-doesnt-kill-triangulations/triangulation_layer0.pdf` — the
   local-survival lemma in two pages.
3. `step-one-kill-statistics/step1_statistics.pdf` — kill rates and
   asymptotics.
