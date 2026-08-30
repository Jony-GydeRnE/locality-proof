#!/usr/bin/env python3
"""Quick check: does any n=10 step-1 survivor need geometric k >= 4?

If not (all die at k=2 or k=3), then K(n) likely saturates at 3.
If any needs k=4, then K(n) keeps growing (structure beyond the pinwheel).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pinwheel_test import (
    all_chords, is_triangulation, kzone_data, survives_geometric_at,
)
from itertools import combinations_with_replacement
import time

n = 10
N = n - 3
print(f"n={n}, multiset size N={N}")

# Step-1 survivors only: skip everything that dies at some k=1 zone
print("scanning for step-1 (k=1) survivors...")
t0 = time.time()
step1_survivors = []
count = 0
for ms in combinations_with_replacement(all_chords(n), N):
    count += 1
    if is_triangulation(ms, n):
        continue
    ms_set = set(ms)
    if all(survives_geometric_at(ms_set, kzone_data(n, 1, i)) for i in range(1, n + 1)):
        step1_survivors.append(ms)
    if count % 200000 == 0:
        print(f"  {count} processed, {len(step1_survivors)} survivors so far ({time.time()-t0:.0f}s)")
print(f"  total enumerated: {count} in {time.time()-t0:.0f}s")
print(f"  step-1 non-tri survivors at n={n}: {len(step1_survivors)}")

# Now find min-kill-k for each
print(f"computing min-kill-k under geometric criterion...")
from collections import Counter
by_k = Counter()
needing_high_k = []
for ms in step1_survivors:
    ms_set = set(ms)
    min_k = None
    for k in range(2, n - 2):
        for i in range(1, n + 1):
            if not survives_geometric_at(ms_set, kzone_data(n, k, i)):
                min_k = k
                break
        if min_k is not None:
            break
    by_k[min_k] += 1
    if min_k is None or min_k >= 4:
        needing_high_k.append((ms, min_k))

print(f"min-kill-k histogram at n={n}:")
for k, c in sorted(by_k.items(), key=lambda x: (x[0] is None, x[0] or 999)):
    print(f"  k = {k}: {c}")

print(f"\nMultisets needing k>=4 or surviving everything: {len(needing_high_k)}")
for ms, k in needing_high_k[:10]:
    print(f"  k={k}: {ms}")
