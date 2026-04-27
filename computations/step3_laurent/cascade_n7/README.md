# Laurent cascade kill at $n = 7$

**Companion to `paper/Laurent series for hard kills/cascade_n7.tex`.**

Symbolic verification (sympy) of the depth-1 Laurent cascade that kills
all seven $n=7$ "fish" multisets — the cyclic family of non-triangulation
multisets that escape Step 1 at every zone.

For each fish $M$, the script:
1. picks a kill zone $\mathcal Z_{r,r+2}$ and Laurent order $k=3$;
2. computes the order-$X_*^{-3}$ coefficient of a chosen free-variable
   fingerprint, producing a 6-term equation in ansatz coefficients;
3. confirms that the five "cousins" of $M$ in that equation are all
   killable by Step 1 at one neighbouring zone;
4. concludes $a_M = 0$.

## Files

- **`cascade_kill_n7.py`** — the verifier. Top-level driver:
  `verify_all_seven_fish()`.
- **`results_cascade_n7.txt`** — full human-readable cascade trace
  (one block per fish: kill zone, order, fingerprint, 6-term equation,
  and the cousins' Step-1 kill zones).

## Run

```bash
python3 cascade_kill_n7.py
```

Runs in seconds. All seven fish kills confirmed.
