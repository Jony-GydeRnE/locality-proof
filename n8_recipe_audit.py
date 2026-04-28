#!/usr/bin/env python3
"""
n8_recipe_audit.py

Parse the 100 n=8 cascade recipes from results_cascade_n8.txt and audit
two structural claims:

  (Claim A) DEPTH-1 UNIFORMITY.  For every survivor, the actual Laurent
            order k matches  ell_kill + 1, where  ell_kill  is the
            number of M-chords that contribute a 1/S factor at the
            chosen kill zone (i.e., specials + bares + non-bare row
            chords).  If true, all 100 fish are depth-1 from leading
            order, irrespective of the absolute k value.

  (Claim B) MINIMUM-ELL ZONE.  The actual kill zone minimizes
            ell across all zones.  If true, the kill-zone selection
            rule is "pick the zone where M has the fewest blocking
            chords."

  (Claim C) COMMON COUSIN ZONE FOR DEPTH-1 (k=3) FISH.  All cousins of
            a depth-1 fish are step-1-killable at one common
            neighbouring zone (the zone whose special equals the
            bare-at-kill of M).  This was the cleanest n=7 feature.

  (Claim D) PHI-LIKE RECIPE STRUCTURE.  fingerprint numerator =
            companion of one of M's non-bare row chords, fingerprint
            denominator = product of M's free chords (with
            multiplicities).
"""

import re
import sys
from collections import Counter

PATH = "/tmp/locality-proof/computations/step3_laurent/cascade_n8/results_cascade_n8.txt"
N_POLY = 8


def normalize(i, j, n=N_POLY):
    i = ((i - 1) % n) + 1
    j = ((j - 1) % n) + 1
    if i > j:
        i, j = j, i
    return (i, j)


def parse_chord(s):
    """Parse '(1,3)' or 'X_1_3' to a tuple."""
    s = s.strip()
    m = re.match(r"\((\d+),(\d+)\)", s)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    m = re.match(r"X_(\d+)_(\d+)", s)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    raise ValueError(f"can't parse chord {s!r}")


def parse_block(block):
    """Extract (M, kill_zone, k, fp_num, fp_denom_chords, cousin_zones) from
    one survivor block of the results file."""
    m_match = re.search(r"M = \{([^}]+)\}", block)
    if not m_match:
        return None
    M_chords = [parse_chord(c) for c in re.findall(r"\(\d+,\d+\)", m_match.group(1))]

    kz_match = re.search(r"Kill zone Z = Z_\{(\d+),(\d+)\}", block)
    kill_zone = (int(kz_match.group(1)), int(kz_match.group(2)))

    k_match = re.search(r"Laurent order k = (\d+)", block)
    k = int(k_match.group(1))

    fp_match = re.search(r"fingerprint = (X_\d+_\d+)/\(([^)]+)\)", block)
    if fp_match is None:
        return None
    fp_num = parse_chord(fp_match.group(1))
    fp_denom_str = fp_match.group(2)
    fp_denom_chords = []
    # Use regex to find each X_a_b possibly followed by **n
    for m in re.finditer(r"X_(\d+)_(\d+)(?:\*\*(\d+))?", fp_denom_str):
        chord = (int(m.group(1)), int(m.group(2)))
        mult = int(m.group(3)) if m.group(3) else 1
        fp_denom_chords.extend([chord] * mult)

    cousin_zones = re.findall(r"step-1 at Z_\{(\d+),(\d+)\}", block)
    cousin_zones = [(int(a), int(b)) for a, b in cousin_zones]

    return {
        "M": M_chords,
        "kill_zone": tuple(sorted(kill_zone)),
        "k": k,
        "fp_num": fp_num,
        "fp_denom": fp_denom_chords,
        "cousin_zones": cousin_zones,
    }


def zone_structure(r, n=N_POLY):
    """At zone Z_{r,r+2} (r is the parameter; note Z is labeled by special chord)."""
    special = normalize(r, r + 2, n)
    bare = normalize(r + 1, r - 1, n)
    pairs = []
    for offset in range(3, n - 1):
        k = ((r - 1 + offset) % n) + 1
        comp = normalize(r, k, n)
        sub = normalize(r + 1, k, n)
        pairs.append((comp, sub))
    return special, bare, pairs


def find_r(zone_special, n=N_POLY):
    """Given special (a,b), find r such that normalize(r, r+2, n) = (a,b)."""
    for r in range(1, n + 1):
        if normalize(r, r + 2, n) == zone_special:
            return r
    return None


def ell_at_zone(M, r, n=N_POLY):
    """Number of M-chords that are bare, non-bare row r+1, or special at zone Z_{r,r+2}."""
    special, bare, pairs = zone_structure(r, n)
    substitutes = [sub for _, sub in pairs]
    blocking = {special, bare} | set(substitutes)
    cnt = 0
    for c in M:
        if c in blocking:
            cnt += 1
    return cnt


def ell_at_kill(M, kill_zone, n=N_POLY):
    r = find_r(kill_zone, n)
    return ell_at_zone(M, r, n)


def step1_killable(M, r, n=N_POLY):
    """K1 + K2: M has no special, no bare, no (companion, substitute) pair."""
    special, bare, pairs = zone_structure(r, n)
    M_set = list(M)
    if M_set.count(special) > 0 or M_set.count(bare) > 0:
        return False
    for comp, sub in pairs:
        if comp in M_set and sub in M_set:
            return False
    return True


# ---- Read all blocks ----
text = open(PATH).read()
# Split at each survivor boundary
sections = re.split(r"(?=Survivor #\d+/)", text)
recipes = []
for s in sections:
    if "Survivor #" in s and "Kill zone" in s:
        rec = parse_block(s)
        if rec is not None:
            recipes.append(rec)

print(f"Parsed {len(recipes)} survivor recipes from {PATH}.")

# ---- CLAIM A: depth-1 from leading? ----
claim_a_match = 0
claim_a_mismatch = []
for rec in recipes:
    ell = ell_at_kill(rec["M"], rec["kill_zone"], N_POLY)
    expected_k = ell + 1
    if rec["k"] == expected_k:
        claim_a_match += 1
    else:
        claim_a_mismatch.append((rec["M"], rec["kill_zone"], rec["k"], ell, expected_k))

print(f"\nCLAIM A (depth-1 from leading: k = ell + 1):")
print(f"  match: {claim_a_match}/100")
if claim_a_mismatch:
    print(f"  mismatches:")
    for M, kz, k, ell, exp_k in claim_a_mismatch[:5]:
        print(f"    M = {M}, kill = {kz}, k = {k}, ell = {ell}, expected k = {exp_k}")

# ---- CLAIM B: kill zone has minimum ell? ----
claim_b_match = 0
claim_b_partial = 0
claim_b_data = []
for rec in recipes:
    M = rec["M"]
    ells = [ell_at_zone(M, r, N_POLY) for r in range(1, N_POLY + 1)]
    # but only over zones where (K1 or K2) violation makes step-1 fail
    relevant_ells = [(r, ell_at_zone(M, r, N_POLY)) for r in range(1, N_POLY + 1)
                     if not step1_killable(M, r, N_POLY)]
    if not relevant_ells:
        claim_b_data.append((M, "NO ZONES BLOCK", -1, -1))
        continue
    min_ell = min(e for _, e in relevant_ells)
    chosen_r = find_r(rec["kill_zone"], N_POLY)
    chosen_ell = ell_at_zone(M, chosen_r, N_POLY)
    if chosen_ell == min_ell:
        claim_b_match += 1
    else:
        claim_b_partial += 1
        claim_b_data.append((M, rec["kill_zone"], chosen_ell, min_ell))

print(f"\nCLAIM B (kill zone minimizes ell over step-1-blocked zones):")
print(f"  match: {claim_b_match}/100   |   non-min ell: {claim_b_partial}/100")
if claim_b_data and claim_b_partial > 0:
    print(f"  example non-min:")
    for M, kz, ce, me in claim_b_data[:3]:
        print(f"    M = {M}, kill = {kz}, chosen ell = {ce}, min ell = {me}")

# ---- CLAIM C: common cousin zone for depth-1 (k=3) recipes? ----
depth1_recipes = [r for r in recipes if r["k"] == 3]
claim_c_match = 0
claim_c_mismatch = 0
for rec in depth1_recipes:
    cousin_z_set = set(rec["cousin_zones"])
    if len(cousin_z_set) == 1:
        claim_c_match += 1
    else:
        claim_c_mismatch += 1

print(f"\nCLAIM C (k=3 fish: all cousins killed at one common zone):")
print(f"  match: {claim_c_match}/{len(depth1_recipes)}   |  multi-zone: {claim_c_mismatch}/{len(depth1_recipes)}")

depth2_recipes = [r for r in recipes if r["k"] == 4]
claim_c2_match = 0
for rec in depth2_recipes:
    cousin_z_set = set(rec["cousin_zones"])
    if len(cousin_z_set) == 1:
        claim_c2_match += 1
print(f"\nFor reference, k=4 fish:")
print(f"  cousins all at one zone: {claim_c2_match}/{len(depth2_recipes)}")
print(f"  (depth-2 cases generally need multiple cousin zones, as expected.)")

# ---- CLAIM D: fingerprint structure ----
# fp_num should be companion of some non-bare row chord of M at kill zone
# fp_denom should be M's free chords (multiset)
claim_d_match = 0
claim_d_mismatch = []
for rec in recipes:
    r = find_r(rec["kill_zone"], N_POLY)
    special, bare, pairs = zone_structure(r, N_POLY)
    substitutes = {sub: comp for comp, sub in pairs}

    M = list(rec["M"])
    # Identify the non-bare row chords in M (i.e., M's chords matching some substitute)
    M_nonbare = [c for c in M if c in substitutes]
    M_companions_avail = [substitutes[c] for c in M_nonbare]
    M_free = [c for c in M if c not in substitutes and c != bare and c != special]

    fp_num_match = rec["fp_num"] in M_companions_avail
    fp_denom_match = sorted(rec["fp_denom"]) == sorted(M_free)
    if fp_num_match and fp_denom_match:
        claim_d_match += 1
    else:
        claim_d_mismatch.append((rec["M"], rec["fp_num"], M_companions_avail,
                                  rec["fp_denom"], M_free, fp_num_match, fp_denom_match))

print(f"\nCLAIM D (fingerprint = X_companion / product(M_free)):")
print(f"  match: {claim_d_match}/100")
if claim_d_mismatch:
    print(f"  mismatches:")
    for M, fpn, comp_avail, fpd, mf, num_ok, denom_ok in claim_d_mismatch[:3]:
        print(f"    M = {M}")
        print(f"      fp_num = {fpn}, available companions = {comp_avail}")
        print(f"      fp_denom = {fpd}, M_free = {mf}")
        print(f"      num_ok={num_ok}, denom_ok={denom_ok}")

# ---- Summary ell distribution ----
ell_dist = Counter(ell_at_kill(r["M"], r["kill_zone"]) for r in recipes)
print(f"\nDistribution of ell_kill (= k - 1):")
for ell, count in sorted(ell_dist.items()):
    print(f"  ell = {ell}  (so k = {ell+1}):  {count} survivors")

# ---- What's the global minimum ell across all step-1-blocked zones? ----
print(f"\n--- min_r ell across step-1-blocked zones (= structural difficulty) ---")
ell_min_dist = Counter()
for rec in recipes:
    M = rec["M"]
    blocked_ells = [ell_at_zone(M, r, N_POLY) for r in range(1, N_POLY + 1)
                    if not step1_killable(M, r, N_POLY)]
    if blocked_ells:
        ell_min_dist[min(blocked_ells)] += 1
for em, count in sorted(ell_min_dist.items()):
    print(f"  min ell = {em}:  {count} survivors")

# ---- Why do some fish need ell_kill > min ell? ----
# Claim B-prime: among zones with min ell, is there ALWAYS one where the
# depth-1 cascade has only step-1-killable cousins?
# (We can't easily check this without re-running the cascade engine; just
#  report the gap.)
print(f"\n--- Zones tried but skipped (chosen ell - min ell > 0) ---")
gap_dist = Counter()
for rec in recipes:
    M = rec["M"]
    blocked_ells = [ell_at_zone(M, r, N_POLY) for r in range(1, N_POLY + 1)
                    if not step1_killable(M, r, N_POLY)]
    if not blocked_ells:
        continue
    min_ell = min(blocked_ells)
    chosen_r = find_r(rec["kill_zone"], N_POLY)
    chosen_ell = ell_at_zone(M, chosen_r, N_POLY)
    gap_dist[chosen_ell - min_ell] += 1
for g, count in sorted(gap_dist.items()):
    print(f"  ell gap = {g}:  {count} survivors")

# ---- Examine the multi-cousin-zone depth-1 cases ----
print(f"\n--- Depth-1 (k=3) fish with cousins killed at multiple zones ---")
mc_count = 0
for rec in depth1_recipes:
    cousin_z_set = set(rec["cousin_zones"])
    if len(cousin_z_set) > 1:
        mc_count += 1
        if mc_count <= 5:
            print(f"  M = {rec['M']}, kill = {rec['kill_zone']},")
            print(f"    cousin zones = {sorted(cousin_z_set)}")
print(f"  ... ({mc_count} total)")

