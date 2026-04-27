"""
================================================================================
analyze_recipes.py  --  Recipe analysis of the n=8 cascade kill results
================================================================================

PURPOSE (BIG PICTURE)
---------------------
The companion script `cascade_kill_n8.py` has verified that all 100 step-1
survivors at n=8 die via a depth-1 Laurent cascade.  This script does the
follow-up structural analysis the n=8 paper needs:

  (1) Parse `results_cascade_n8.txt` to extract per-survivor (M, kill zone,
      Laurent order, companion).

  (2) Compute the cyclic-orbit decomposition of the 100 survivors under the
      Z_8 action  v -> v + s (mod 8).  Each orbit has size dividing 8.

  (3) For each survivor, compute the missed-vertex set
      V_missed(M) = {1..8} \ V(M),  where V(M) = vertices used by M.
      List the cyclic-distance-2 pairs inside V_missed (the "frame"
      candidates -- specials that could be picked as the kill zone).

  (4) Compare the empirical kill zone Z_{r, r+2} from the cascade run with
      the Phi-I prediction: a kill zone whose special chord (r, r+2) has
      both endpoints in V_missed (i.e., the frame is supported by missed
      vertices).  Report the Phi-I match rate.

OUTPUTS
-------
- recipe_analysis.md    -- markdown table per survivor (M, V_missed,
                            frame candidates, actual kill zone, Phi-I match,
                            orbit ID).
- orbits.md             -- summary: number of cyclic orbits, sizes, and one
                            canonical representative per orbit.
- analysis_summary.txt  -- short end-of-run summary (Phi-I match rate,
                            orbit count, etc.).

CODE STYLE
----------
Per repo contribution standards (see root README), every function has a
docstring with both a LOGIC section (what the algorithm does in CS terms)
and a PHYSICS / MATHEMATICS section (what it computes in proof language).
"""

from collections import defaultdict
import re
import os


# ============================================================================
# 1.  Parse `results_cascade_n8.txt`
# ============================================================================

RESULTS_FILE = os.path.join(os.path.dirname(__file__), "results_cascade_n8.txt")


def parse_results(path=RESULTS_FILE):
    """
    Parse the cascade-run results file into a list of per-survivor dicts.

    LOGIC:
      Scan the file line by line for:
        "Survivor #N/100: M = {(i,j),...}"
        "  Kill zone Z = Z_{r,r+2};  Laurent order k = K;"
        "  fingerprint = ... (companion (i,j))."
      Group these into a single record per survivor.

    PHYSICS:  the records are exactly the (multiset, kill recipe) pairs
              that the cascade verifier produced.  Each record represents
              a step-1 survivor M and its depth-1 Laurent cascade kill at
              zone Z_{r,r+2} at order X_*^{-k} with fingerprint involving
              the listed companion chord.

    Returns a list of dicts with keys:
        index, multiset, kill_r, kill_order, companion.
    """
    survivors = []
    current = None
    with open(path, "r") as f:
        for line in f:
            m_idx = re.match(
                r"^Survivor #(\d+)/\d+: M = \{(.+)\}\s*$", line)
            if m_idx:
                if current is not None:
                    survivors.append(current)
                idx = int(m_idx.group(1))
                multiset = tuple(
                    tuple(int(x) for x in pair.strip("()").split(","))
                    for pair in re.findall(r"\(\d+,\d+\)", m_idx.group(2))
                )
                current = {"index": idx, "multiset": multiset,
                           "kill_r": None, "kill_order": None,
                           "companion": None}
                continue
            m_kill = re.match(
                r"^\s+Kill zone Z = Z_\{(\d+),(\d+)\};\s+"
                r"Laurent order k = (\d+);\s*$", line)
            if m_kill and current is not None:
                # The "kill_r" we record is the smaller index of the pair,
                # i.e., the r in zone Z_{r, r+2}.  The pair (r, r+2) on the
                # 8-gon is normalized so r < r+2 mod 8.
                a, b = int(m_kill.group(1)), int(m_kill.group(2))
                # Identify which is the "r" (the eliminated-vertex zone
                # label).  Convention: zone_r = a iff (a, a+2) == (a, b)
                # cyclically.  At n=8, one of (a, b-a) is (r, 2) cyclically.
                # We just store both and recover (r, r+2) downstream.
                current["kill_zone"] = (a, b)
                current["kill_order"] = int(m_kill.group(3))
                continue
            m_fp = re.match(
                r"^\s+fingerprint = .* \(companion \((\d+), (\d+)\)\)\.\s*$",
                line)
            if m_fp and current is not None:
                current["companion"] = (int(m_fp.group(1)),
                                        int(m_fp.group(2)))
                continue
        if current is not None:
            survivors.append(current)
    return survivors


# ============================================================================
# 2.  Vertex / missed-vertex structure
# ============================================================================

def vertices_of(multiset):
    """
    The set of polygon vertices that appear as endpoints of any chord in M.

    LOGIC:  flatten all chord pairs into a set of integers.

    PHYSICS:  V(M) records which vertices the multiset's denominator
              "touches"; a vertex r+1 not in V(M) means M cannot contain
              the bare (r+1, r-1) or any non-bare substitute X_{r+1,k} on
              zone Z_{r,r+2}, simplifying the killability analysis.
    """
    vs = set()
    for (i, j) in multiset:
        vs.add(i)
        vs.add(j)
    return vs


def missed_vertices(multiset, n=8):
    """
    V_missed(M) = {1, ..., n} \\ V(M).

    LOGIC:  set difference.

    PHYSICS:  vertices not touched by any chord of M.  These are the
              candidate eliminated-vertex labels r+1 for zones at which
              M trivially satisfies (K1) and (K2) for all (companion,
              substitute) pairs at that vertex -- so M is closer to being
              step-1 killable at those zones (though the special and bare
              chords may still block it).
    """
    return set(range(1, n + 1)) - vertices_of(multiset)


def cyclic_distance(a, b, n=8):
    """
    Minimum cyclic distance between vertices a, b on the n-gon.

    LOGIC:  d = min(|a - b|, n - |a - b|) in {0, 1, ..., n//2}.

    PHYSICS:  on the 8-gon, distance-2 pairs (a, a+2 mod 8) are exactly
              the pairs that are endpoints of a "special" chord X_{a,a+2}
              at some zone Z_{a, a+2}; these are the candidates for the
              Phi-I-predicted kill zone.
    """
    d = abs(a - b) % n
    return min(d, n - d)


def frame_pairs(missed, n=8):
    """
    All pairs (a, b) inside V_missed at cyclic distance 2.

    LOGIC:  iterate over ordered pairs in missed; keep those with cyclic
            distance exactly 2.

    PHYSICS:  each such pair (a, a+2 mod n) corresponds to the special
              chord X_{a, a+2} of zone Z_{a, a+2}; if this special chord
              is not in M and the bare X_{a+1, a-1} is not in M, then
              (K1) holds at zone Z_{a, a+2} from the "frame" structure
              alone.  So these are the candidate Phi-I kill zones.
    """
    out = []
    missed = sorted(missed)
    for i, a in enumerate(missed):
        for b in missed[i + 1:]:
            if cyclic_distance(a, b, n) == 2:
                # Canonical zone label (r, r+2): pick the one with r+2
                # cyclically after r.  Both (a, b) and (b, a) give the
                # same zone, normalised by smaller r modulo n.
                # The chord (a, b) on the n-gon is normalized with a < b,
                # but the zone label is (r, r+2) where (r+2) - r = 2 mod n.
                # So if b - a == 2, zone = (a, b); if a - b == -2 mod n
                # (i.e. a > b and n - (a - b) == 2), zone wraps cyclically.
                if (b - a) % n == 2:
                    out.append((a, b))
                else:
                    out.append((b, a))  # b is the "r", wraps to a
    return out


# ============================================================================
# 3.  Cyclic orbit decomposition
# ============================================================================

def shift_chord(c, s, n=8):
    """
    Apply the cyclic shift v -> v + s (mod n) to a chord c = (i, j).

    LOGIC:  add s to each endpoint, reduce mod n into {1..n}, sort.

    PHYSICS:  the Z_n action on chord multisets.  Two multisets are in
              the same orbit iff they differ by a cyclic shift.
    """
    i, j = c
    i2 = ((i - 1 + s) % n) + 1
    j2 = ((j - 1 + s) % n) + 1
    if i2 > j2:
        i2, j2 = j2, i2
    return (i2, j2)


def shift_multiset(M, s, n=8):
    """
    Apply v -> v + s to every chord of M, return the sorted-tuple form.

    LOGIC:  map shift_chord across M, then sort to get a canonical hash.

    PHYSICS:  this is the action of Z_n on the chord-multiset ansatz;
              two multisets are cyclically equivalent iff some shift maps
              one to the other.
    """
    return tuple(sorted(shift_chord(c, s, n) for c in M))


def canonical_orbit_rep(M, n=8):
    """
    Lexicographically smallest cyclic-shift image of M.

    LOGIC:  compute shift_multiset(M, s) for s in 0..n-1, return the min.

    PHYSICS:  the chosen orbit representative (one per Z_n-orbit).  Two
              multisets are in the same orbit iff their canonical reps
              are identical.
    """
    return min(shift_multiset(M, s, n) for s in range(n))


def orbit_decomposition(survivors, n=8):
    """
    Partition the survivors into Z_n-orbits.

    LOGIC:  group survivors by canonical_orbit_rep.  Orbit sizes divide n.

    PHYSICS:  the cyclic action on the ansatz preserves the kill mechanism
              (cyclic shifts of zones produce cyclic shifts of recipes), so
              two survivors in the same orbit have related kill recipes
              and one orbit-representative analysis gives the kill recipe
              for the whole orbit.
    """
    orbits = defaultdict(list)
    for s in survivors:
        rep = canonical_orbit_rep(s["multiset"], n)
        orbits[rep].append(s)
    return orbits


# ============================================================================
# 4.  Phi-I prediction
# ============================================================================

def phi_i_predicted_zones(multiset, n=8):
    """
    The set of kill zones predicted by the Phi-I rule from the missed-vertex
    structure of M.

    LOGIC:  return frame_pairs(missed_vertices(M)).  Each (r, r+2 cyclic)
            in the result is a predicted kill zone Z_{r, r+2}.

    PHYSICS:  Phi-I rule: M is killable by a depth-1 Laurent cascade at
              zone Z_{r, r+2} whenever the special chord (r, r+2) has both
              endpoints in V_missed(M) -- i.e., the special is not in M
              and the bare X_{r+1, r-1} is also not in M (because vertex
              r+1 is in V_missed too -- though here we only require
              endpoints of the special itself to be missed; the full
              kill recipe requires more).  This is the "frame
              configuration" hypothesis: the frame {r, r+2} sits inside
              V_missed.
    """
    missed = missed_vertices(multiset, n)
    return frame_pairs(missed, n)


def phi_i_match(survivor, n=8):
    """
    Does the empirical kill zone of `survivor` match a Phi-I prediction?

    LOGIC:  compute predicted zones; check membership of the actual
            kill_zone.

    PHYSICS:  if Phi-I is correct as a structural rule, every survivor's
              empirical kill zone should be in its predicted-zone list.
    """
    pred = phi_i_predicted_zones(survivor["multiset"], n)
    return survivor["kill_zone"] in pred


# ============================================================================
# 5.  Reporting
# ============================================================================

def format_multiset(M):
    """Pretty-print a multiset as {(i,j),...}."""
    return "{" + ",".join(f"({i},{j})" for (i, j) in M) + "}"


def write_recipe_analysis(survivors, orbits, output_path,
                          analysis_summary_path):
    """
    Write the per-survivor analysis table + summary.

    LOGIC:  one markdown table row per survivor with columns: index,
            multiset, missed vertices, frame candidates, actual kill
            zone, Phi-I match, orbit ID.  Plus an end-of-file summary.

    PHYSICS:  the table is the "raw analytical data" the n=8 paper will
              cite.  Each row tells you for one survivor whether Phi-I's
              missed-vertex prediction matches the empirical recipe.
    """
    # Number orbits 1..N by sorted canonical rep order.
    orbit_reps_sorted = sorted(orbits.keys())
    orbit_id_map = {rep: i + 1 for i, rep in enumerate(orbit_reps_sorted)}

    n = 8
    matches = 0
    lines = []
    lines.append("# n=8 Step-1-survivor recipe analysis\n")
    lines.append(f"All {len(survivors)} step-1 survivors at $n=8$ are killed "
                 f"by a depth-1 Laurent cascade (verified by "
                 f"`cascade_kill_n8.py`). This file analyses the structural "
                 f"recipe for each one.\n")
    lines.append("## Per-survivor table\n")
    lines.append("| # | M | V_missed | Frame candidates "
                 "(Phi-I predicted zones) | Actual zone | Order | "
                 "Phi-I match? | Orbit |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for s in survivors:
        M = s["multiset"]
        miss = sorted(missed_vertices(M, n))
        pred = phi_i_predicted_zones(M, n)
        actual = s["kill_zone"]
        match = (actual in pred)
        if match:
            matches += 1
        rep = canonical_orbit_rep(M, n)
        oid = orbit_id_map[rep]
        miss_str = "{" + ",".join(str(v) for v in miss) + "}"
        pred_str = ",".join(f"({a},{b})" for (a, b) in pred) or "(none)"
        match_str = "YES" if match else "no"
        lines.append(f"| {s['index']} | {format_multiset(M)} | {miss_str} | "
                     f"{pred_str} | {actual} | {s['kill_order']} | "
                     f"{match_str} | {oid} |")

    lines.append("")
    lines.append("## Phi-I match rate\n")
    lines.append(f"- {matches} / {len(survivors)} = "
                 f"{100.0 * matches / len(survivors):.1f}%  --- empirical "
                 f"kill zone matches a Phi-I-predicted (frame-configuration) "
                 f"zone.\n")

    lines.append("## Orbit summary\n")
    sizes = [len(v) for v in orbits.values()]
    size_dist = defaultdict(int)
    for s in sizes:
        size_dist[s] += 1
    lines.append(f"- Total orbits: {len(orbits)}.")
    lines.append(f"- Sum of orbit sizes: {sum(sizes)} (expected: "
                 f"{len(survivors)}).\n")
    lines.append("Orbit-size distribution:\n")
    lines.append("| Orbit size | # orbits |")
    lines.append("|---|---:|")
    for sz in sorted(size_dist.keys()):
        lines.append(f"| {sz} | {size_dist[sz]} |")
    lines.append("")
    lines.append("Per-orbit summary (canonical representative, size, "
                 "Phi-I match count within orbit):\n")
    lines.append("| Orbit | Canonical rep | Size | Phi-I match in orbit |")
    lines.append("|---|---|---:|---|")
    for rep in orbit_reps_sorted:
        oid = orbit_id_map[rep]
        members = orbits[rep]
        in_orbit_match = sum(1 for s in members if phi_i_match(s, n))
        lines.append(f"| {oid} | {format_multiset(rep)} | "
                     f"{len(members)} | {in_orbit_match} / {len(members)} |")

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    # Also write a short end-summary plain-text file.
    with open(analysis_summary_path, "w") as f:
        f.write(
            f"n=8 step-1-survivor recipe analysis -- summary\n"
            f"=============================================\n"
            f"\n"
            f"Total survivors:          {len(survivors)}\n"
            f"Total cyclic orbits:      {len(orbits)}\n"
            f"Phi-I match rate:         {matches}/{len(survivors)} "
            f"= {100.0 * matches / len(survivors):.1f}%\n"
            f"\n"
            f"Phi-I rule:  the empirical kill zone Z_{{r,r+2}} of M\n"
            f"             satisfies (r, r+2) in V_missed(M), i.e., both\n"
            f"             endpoints of the special chord are vertices\n"
            f"             that no chord of M touches.\n"
            f"\n"
            f"Orbit-size distribution:\n"
        )
        for sz in sorted(size_dist.keys()):
            f.write(f"  size {sz}: {size_dist[sz]} orbit(s)\n")
        f.write("\nFull per-survivor table: see recipe_analysis.md\n")


def write_orbits_md(orbits, output_path, n=8):
    """
    Write the orbit-decomposition summary.

    LOGIC:  list each orbit, its canonical representative, its size, and
            the index of one of its members (so the reader can find the
            full record in the cascade trace).

    PHYSICS:  the orbit decomposition is the right structural object for
              the all-n proof: the kill recipe for one orbit member
              determines the kill recipe for the whole orbit by cyclic
              shift.  So in principle one needs to verify the recipe
              once per orbit, not once per survivor.
    """
    orbit_reps_sorted = sorted(orbits.keys())
    lines = []
    lines.append(f"# n={n} step-1-survivor cyclic orbits\n")
    lines.append(f"Group: Z_{n} acting by  v -> v + s (mod {n}).  "
                 f"Orbit sizes divide {n}.\n")
    lines.append(f"Total orbits: {len(orbits)}.")
    sizes = sorted({len(v) for v in orbits.values()})
    lines.append(f"Distinct orbit sizes: {sizes}\n")
    lines.append("| Orbit | Size | Canonical rep | One member index |")
    lines.append("|---:|---:|---|---:|")
    for i, rep in enumerate(orbit_reps_sorted):
        members = orbits[rep]
        one_idx = min(s["index"] for s in members)
        lines.append(f"| {i+1} | {len(members)} | {format_multiset(rep)} | "
                     f"{one_idx} |")
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ============================================================================
# 6.  Main
# ============================================================================

def main():
    """
    Top-level driver: parse, analyse, write.

    LOGIC:  parse_results -> orbit_decomposition -> write all three output
            files (recipe_analysis.md, orbits.md, analysis_summary.txt).

    PHYSICS:  produces the "raw analytical data" that the n=8 follow-up
              paper will cite.  All three outputs are human-readable.
    """
    survivors = parse_results()
    n = 8
    if not survivors:
        print("ERROR: failed to parse any survivors from results file.")
        return
    print(f"Parsed {len(survivors)} survivors from {RESULTS_FILE}.")

    orbits = orbit_decomposition(survivors, n)
    print(f"Cyclic-orbit decomposition: {len(orbits)} orbits, sizes "
          f"{sorted({len(v) for v in orbits.values()})}.")

    here = os.path.dirname(__file__)
    write_recipe_analysis(
        survivors, orbits,
        output_path=os.path.join(here, "recipe_analysis.md"),
        analysis_summary_path=os.path.join(here, "analysis_summary.txt"),
    )
    write_orbits_md(orbits, output_path=os.path.join(here, "orbits.md"),
                    n=n)

    matches = sum(1 for s in survivors if phi_i_match(s, n))
    print(f"Phi-I match rate: {matches}/{len(survivors)} = "
          f"{100.0 * matches / len(survivors):.1f}%")
    print(f"Outputs: recipe_analysis.md, orbits.md, analysis_summary.txt")


if __name__ == "__main__":
    main()
