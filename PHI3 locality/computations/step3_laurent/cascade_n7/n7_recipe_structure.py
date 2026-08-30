#!/usr/bin/env python3
"""
n7_recipe_structure.py

Take the seven n=7 fish recipes from cascade_kill_n7.py and decompose each
into structural primitives:

  - frame vertex  v*    (the shared vertex of the two framing chords)
  - framing chords      (the two M-chords incident to v*)
  - crossing pair       (the two M-chords that cross; both avoid v*)
  - kill zone           Z_{v*-1, v*+1}   (eliminates v*)
  - bare-at-kill        the framing chord {v*, v*-2}
  - non-bare-at-kill    the framing chord {v*, v*+2}     (well, mod n)
  - companion           X_{v*-1, k_M}    (where k_M = non-shared endpoint
                                          of the non-bare framing chord)
  - fingerprint         X_{v*-1, k_M} / (X_a X_b)
                        where (X_a, X_b) is the crossing pair
  - cousin zone         Z_{v*-2, v*}     (special = bare-at-kill)

Then verify that the recipe predicted by Phi(M) matches the recipe in
results_cascade_n7.txt.

This is the cleanest characterization of the depth-1 cascade map at n=7,
suitable for testing against n=8 data once the n=8 cascade run completes.
"""

from itertools import combinations

N_POLY = 7  # heptagon


def normalize(i, j, n=N_POLY):
    """Cyclic-canonical chord name (a, b) with 1 <= a < b <= n."""
    i = ((i - 1) % n) + 1
    j = ((j - 1) % n) + 1
    if i > j:
        i, j = j, i
    return (i, j)


def chords_cross(c1, c2, n=N_POLY):
    """Two chords cross iff their endpoints alternate on the cycle."""
    a, b = c1
    c, d = c2
    if {a, b} & {c, d}:
        return False
    inside = lambda v: a < v < b
    return inside(c) != inside(d)


def find_crossings(M, n=N_POLY):
    """Return list of crossing pairs (c1, c2) in multiset M."""
    out = []
    for i, j in combinations(range(len(M)), 2):
        if chords_cross(M[i], M[j], n):
            out.append((M[i], M[j]))
    return out


def find_frame_vertex(M, n=N_POLY):
    """
    Find vertex v* incident to (at least) two chords of M, such that the
    OTHER two chords are a crossing pair avoiding v*.

    Returns (v_star, framing_chords, crossing_pair) or None if no such v*.
    """
    for v in range(1, n + 1):
        framing = [c for c in M if v in c]
        non_framing = [c for c in M if v not in c]
        if len(framing) != 2 or len(non_framing) != 2:
            continue
        if chords_cross(non_framing[0], non_framing[1], n):
            return v, tuple(framing), tuple(non_framing)
    return None


def predict_recipe(M, n=N_POLY):
    """Apply the conjectured recipe Phi to multiset M; return a dict."""
    info = find_frame_vertex(M, n)
    if info is None:
        return None
    v_star, framing, crossing = info

    # Identify bare and non-bare framing chords.
    # bare framing chord = (v*-2, v*),  non-bare = (v*, v*+2)
    bare_chord = normalize(v_star - 2, v_star, n)
    nonbare_chord = normalize(v_star, v_star + 2, n)

    # but framing chords may not match these exactly — find them in M.
    if bare_chord not in framing or nonbare_chord not in framing:
        return None  # frame vertex doesn't have the canonical (v*-2, v*+2) framing

    # k_M = non-shared endpoint of nonbare framing chord
    k_M = nonbare_chord[0] if nonbare_chord[0] != v_star else nonbare_chord[1]

    # companion = X_{v*-1, k_M}
    r = (v_star - 1 - 1) % n + 1  # r such that r+1 = v*  (so r = v*-1)
    companion = normalize(r, k_M, n)

    # kill zone = Z_{r, r+2} with r = v*-1
    kill_zone = normalize(r, r + 2, n)

    # cousin zone = Z_{v*-2, v*}  (special = bare-at-kill)
    cousin_zone = bare_chord

    return {
        "frame_vertex": v_star,
        "framing_chords": framing,
        "crossing_pair": crossing,
        "bare_chord": bare_chord,
        "nonbare_chord": nonbare_chord,
        "k_M": k_M,
        "companion": companion,
        "kill_zone": kill_zone,
        "cousin_zone": cousin_zone,
        "fingerprint_num": companion,
        "fingerprint_denom": crossing,
    }


# ----- The seven n=7 fish from results_cascade_n7.txt -----
fish_data = [
    ("Fish #1", ((1, 3), (1, 6), (2, 4), (4, 6)),
                {"kill_zone": (5, 7), "fingerprint_num": (1, 5),
                 "fingerprint_denom": ((1, 3), (2, 4)), "cousin_zone": (4, 6)}),
    ("Fish #2", ((1, 3), (1, 6), (3, 5), (4, 6)),
                {"kill_zone": (2, 7), "fingerprint_num": (3, 7),
                 "fingerprint_denom": ((3, 5), (4, 6)), "cousin_zone": (1, 6)}),
    ("Fish #3", ((1, 3), (1, 6), (3, 5), (5, 7)),
                {"kill_zone": (2, 4), "fingerprint_num": (2, 5),
                 "fingerprint_denom": ((1, 6), (5, 7)), "cousin_zone": (1, 3)}),
    ("Fish #4", ((1, 3), (2, 7), (3, 5), (5, 7)),
                {"kill_zone": (4, 6), "fingerprint_num": (4, 7),
                 "fingerprint_denom": ((1, 3), (2, 7)), "cousin_zone": (3, 5)}),
    ("Fish #5", ((1, 6), (2, 4), (2, 7), (4, 6)),
                {"kill_zone": (3, 5), "fingerprint_num": (3, 6),
                 "fingerprint_denom": ((1, 6), (2, 7)), "cousin_zone": (2, 4)}),
    ("Fish #6", ((2, 4), (2, 7), (3, 5), (5, 7)),
                {"kill_zone": (1, 6), "fingerprint_num": (2, 6),
                 "fingerprint_denom": ((2, 4), (3, 5)), "cousin_zone": (5, 7)}),
    ("Fish #7", ((2, 4), (2, 7), (4, 6), (5, 7)),
                {"kill_zone": (1, 3), "fingerprint_num": (1, 4),
                 "fingerprint_denom": ((4, 6), (5, 7)), "cousin_zone": (2, 7)}),
]


def fmt_chord(c):
    return f"({c[0]},{c[1]})"


def fmt_pair(p):
    a, b = sorted(p)
    return f"{{{fmt_chord(a)},{fmt_chord(b)}}}"


print("=" * 78)
print("Verifying Phi(M) recipe predictions against cascade_n7 actual results")
print("=" * 78)

all_match = True
for name, M, actual in fish_data:
    pred = predict_recipe(M, N_POLY)
    print(f"\n{name}:  M = {{ {', '.join(fmt_chord(c) for c in M)} }}")
    if pred is None:
        print("  Phi failed: no frame vertex with crossing pair")
        all_match = False
        continue

    print(f"  frame vertex v* = {pred['frame_vertex']}")
    print(f"  framing chords  = {fmt_pair(pred['framing_chords'])}")
    print(f"  crossing pair   = {fmt_pair(pred['crossing_pair'])}")
    print(f"  bare-at-kill    = {fmt_chord(pred['bare_chord'])}")
    print(f"  non-bare        = {fmt_chord(pred['nonbare_chord'])}")
    print(f"  k_M             = {pred['k_M']}")
    print(f"  companion       = {fmt_chord(pred['companion'])}")

    actual_denom = tuple(sorted(actual["fingerprint_denom"]))
    pred_denom = tuple(sorted(pred["fingerprint_denom"]))
    matches = (
        pred["kill_zone"] == actual["kill_zone"]
        and pred["fingerprint_num"] == actual["fingerprint_num"]
        and actual_denom == pred_denom
        and pred["cousin_zone"] == actual["cousin_zone"]
    )
    flag = "MATCH" if matches else "MISMATCH"
    print(f"  predicted: kill={fmt_chord(pred['kill_zone'])}  "
          f"fp_num={fmt_chord(pred['fingerprint_num'])}  "
          f"fp_denom={fmt_pair(pred['fingerprint_denom'])}  "
          f"cousin={fmt_chord(pred['cousin_zone'])}")
    print(f"  actual:    kill={fmt_chord(actual['kill_zone'])}  "
          f"fp_num={fmt_chord(actual['fingerprint_num'])}  "
          f"fp_denom={fmt_pair(actual['fingerprint_denom'])}  "
          f"cousin={fmt_chord(actual['cousin_zone'])}")
    print(f"  -> {flag}")
    if not matches:
        all_match = False

print()
print("=" * 78)
if all_match:
    print("VERIFIED: Phi correctly predicts all 7 n=7 fish recipes from local")
    print("structural primitives (frame vertex + framing chords + crossing pair).")
else:
    print("Phi prediction failed on at least one fish.")
print("=" * 78)
