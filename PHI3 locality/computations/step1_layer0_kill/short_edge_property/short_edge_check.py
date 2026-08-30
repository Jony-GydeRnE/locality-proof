"""
Conjecture check: do all step-1 non-tri survivors have a "short, isolated" chord?
================================================================================

THE QUESTION
------------
A chord {a, b} on the n-gon has LENGTH = min(|a-b|, n-|a-b|), so length ∈
{2, 3, ..., floor(n/2)}.  Length 1 = polygon edge (not a chord); length 2
is the SHORTEST possible chord.  Two chords CROSS iff their endpoints
alternate around the polygon boundary.

We say a multiset M has a "short isolated chord" if there exists a chord
c ∈ M with length 2 such that c does NOT cross any other chord in M.
(Sharing an endpoint is fine and doesn't count as crossing.)

This script asks, for each n in {7, 8, 9}:

      Does EVERY step-1 non-tri survivor have a short isolated chord?

If yes for n=8 and n=9, that gives a free, easy-to-state structural
property of the surviving multisets — every survivor has at least one
length-2 chord that "lives alone" in the diagram.

WHAT IT DOES
------------
1. Enumerate all size-(n-3) multisets of chords on the n-gon.
2. Filter to step-1 survivors using the standard zone-by-zone kill rule
   (re-implemented here to keep this script self-contained).
3. Drop triangulations (we only care about the non-tri survivors that
   feed the step-2/step-3 kill machinery).
4. For each survivor, check the short-isolated-chord property.
5. Report counts and counterexamples (if any).

OUTPUTS
-------
  computations/step1_layer0_kill/short_edge_property/
      short_edge_check.json   -- per-n results + counterexample list
      short_edge_check.md     -- human-readable summary
"""

import json
import os
import sys
from itertools import combinations_with_replacement


# ============================================================================
# 1.  CHORDS, ZONES, AND THE STEP-1 KILL CRITERION (self-contained)
# ============================================================================
# These are the same definitions as step1_uncaught.py; copied here so this
# script doesn't have to import across the nested folder structure.

def normalize(i, j, n):
    """Canonical (a, b) with a < b in {1..n}, allowing cyclic input."""
    i = ((i - 1) % n) + 1
    j = ((j - 1) % n) + 1
    if i > j:
        i, j = j, i
    return (i, j)


def all_chords(n):
    """All planar chords (= non-edge unordered vertex pairs) of the n-gon."""
    return [(i, j) for i in range(1, n + 1)
                   for j in range(i + 1, n + 1)
                   if j - i >= 2 and not (i == 1 and j == n)]


def zone_structure(r, n):
    """At zone Z_{r,r+2}: special chord, bare chord, and (companion, substitute)
       pairs.  See step1_uncaught.py for the derivation."""
    special = normalize(r, r + 2, n)
    bare    = normalize(r + 1, r - 1, n)
    pairs = []
    for offset in range(3, n - 1):
        k = ((r - 1 + offset) % n) + 1
        companion  = normalize(r,     k, n)
        substitute = normalize(r + 1, k, n)
        pairs.append((companion, substitute))
    return special, bare, pairs


def killable_at_zone(ms_set, structure):
    """The step-1 kill criterion at a single zone.  ms_set is a Python set
       (NOT multiset) of chords; ms_set ⊂ all_chords(n).  Returns True iff
       this zone kills the multiset's coefficient under the leading-order
       Laurent argument."""
    special, bare, pairs = structure
    if special in ms_set or bare in ms_set:
        return False
    for companion, substitute in pairs:
        if companion in ms_set and substitute in ms_set:
            return False
    return True


def chords_cross(c1, c2):
    """Two chords cross iff their endpoints alternate around the polygon
       boundary (and they don't share an endpoint)."""
    a, b = c1
    c, d = c2
    if {a, b} & {c, d}:
        return False
    inside = lambda v: a < v < b
    return inside(c) != inside(d)


def is_triangulation(ms):
    """Distinct chords with no two of them crossing — these are the
       Catalan-many tree-amplitude diagrams."""
    if len(set(ms)) != len(ms):
        return False
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            if chords_cross(ms[i], ms[j]):
                return False
    return True


# ============================================================================
# 2.  CHORD LENGTH AND THE "SHORT ISOLATED" PROPERTY
# ============================================================================

def chord_length(c, n):
    """Cyclic distance between the two endpoints; 2 = shortest chord."""
    a, b = c
    d = abs(b - a)
    return min(d, n - d)


def has_short_isolated_chord(multiset, n):
    """Return (True, witness_chord) if some length-2 chord c ∈ multiset
       does NOT cross any other chord in the multiset; otherwise (False, None).
       Sharing an endpoint with another chord is fine — that's not a crossing."""
    for i, c in enumerate(multiset):
        if chord_length(c, n) != 2:
            continue
        # Test c against every other chord (including duplicate copies of
        # the same chord-pair, which by definition share endpoints and
        # therefore don't cross).
        crosses_something = False
        for j, other in enumerate(multiset):
            if j == i:
                continue
            if chords_cross(c, other):
                crosses_something = True
                break
        if not crosses_something:
            return True, c
    return False, None


# ============================================================================
# 3.  ENUMERATE STEP-1 NON-TRI SURVIVORS AND CHECK THE PROPERTY
# ============================================================================

def survivors_and_check(n, verbose=True):
    chords = all_chords(n)
    structures = [zone_structure(r, n) for r in range(1, n + 1)]
    N = n - 3   # multiset size
    if verbose:
        print(f"\n--- n = {n} : enumerating size-{N} multisets of "
              f"{len(chords)} chords ---")

    total = 0
    n_uncaught = 0
    n_nontri = 0
    n_with = 0          # # nontri-survivors with a short isolated chord
    n_without = 0       # # nontri-survivors without one
    counterexamples = []

    for ms in combinations_with_replacement(chords, N):
        total += 1
        ms_set = set(ms)
        # step-1: must be uncaught at EVERY zone
        if any(killable_at_zone(ms_set, s) for s in structures):
            continue
        n_uncaught += 1
        # only non-triangulation survivors are interesting
        if is_triangulation(list(ms)):
            continue
        n_nontri += 1
        ok, witness = has_short_isolated_chord(list(ms), n)
        if ok:
            n_with += 1
        else:
            n_without += 1
            if len(counterexamples) < 50:
                counterexamples.append([list(c) for c in ms])

    if verbose:
        print(f"  total multisets    : {total}")
        print(f"  step-1 uncaught    : {n_uncaught}")
        print(f"  non-tri survivors  : {n_nontri}")
        print(f"  WITH short-isol-chord    : {n_with}")
        print(f"  WITHOUT (counterexamples): {n_without}")
        if n_without:
            print(f"  showing first {min(5, len(counterexamples))} counterexamples:")
            for ce in counterexamples[:5]:
                print(f"    {ce}")

    return {
        "n": n,
        "total_multisets": total,
        "step1_uncaught": n_uncaught,
        "nontri_survivors": n_nontri,
        "with_short_isolated_chord": n_with,
        "without_short_isolated_chord": n_without,
        "all_have_property": (n_without == 0),
        "counterexamples": counterexamples,
    }


# ============================================================================
# 4.  MAIN
# ============================================================================

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "short_edge_check.json")
OUT_MD   = os.path.join(HERE, "short_edge_check.md")


def main():
    results = []
    for n in (7, 8, 9):
        results.append(survivors_and_check(n))

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote: {OUT_JSON}")

    # Markdown summary
    lines = []
    lines.append("# Does every step-1 non-tri survivor have a short isolated chord?\n")
    lines.append("**Property checked:** ∃ chord c ∈ multiset M with length 2 such "
                 "that c does NOT cross any other chord in M.  (Length 2 is the "
                 "shortest chord on an n-gon.  Sharing an endpoint with another "
                 "chord does not count as crossing.)\n")
    lines.append("## Per-n results\n")
    lines.append("| n | non-tri survivors | with property | without | universal? |")
    lines.append("|---:|---:|---:|---:|:---:|")
    for r in results:
        verdict = "✅ YES" if r["all_have_property"] else "❌ NO"
        lines.append(f"| {r['n']} | {r['nontri_survivors']} | "
                     f"{r['with_short_isolated_chord']} | "
                     f"{r['without_short_isolated_chord']} | {verdict} |")
    lines.append("")

    for r in results:
        if not r["all_have_property"]:
            lines.append(f"## n = {r['n']} — counterexamples\n")
            lines.append(f"{r['without_short_isolated_chord']} multisets have NO "
                         f"length-2 chord that avoids crossings.  First "
                         f"{min(20, len(r['counterexamples']))} shown:\n")
            for ce in r["counterexamples"][:20]:
                ms_str = '{' + ', '.join(f"({c[0]},{c[1]})" for c in ce) + '}'
                lines.append(f"- `{ms_str}`")
            lines.append("")

    with open(OUT_MD, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote: {OUT_MD}")


if __name__ == "__main__":
    main()
