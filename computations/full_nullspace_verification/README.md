# Full-ansatz nullspace verification

Symbolic / numerical computation that the full constraint system

  (most general rational ansatz B with simple planar poles) ∩ (B vanishes on every cyclic 1-zero locus)

has **nullspace dimension exactly 1**, and the surviving direction equals A_n^tree.

This is the strongest possible verification of the locality + unitarity theorem at finite n.

## Files
- `full_ansatz_n5.py` — n=5 via SymPy (15-dim ansatz, 5 loci, ~10 s).
- `full_ansatz_n6.py` — n=6 numerical (84- or 165-dim ansatz, < 5 s).
- `full_ansatz_n7.py` — n=7 numerical.
- `results_n5.txt`, `results_n6.txt`, `results_n7.txt` — saved console outputs.

## Headline result
At every n ∈ {5, 6, 7}, the nullspace is one-dimensional and is a scalar multiple of A_n^tree.

## Run
```bash
python3 full_ansatz_n5.py
python3 full_ansatz_n6.py
python3 full_ansatz_n7.py
```
