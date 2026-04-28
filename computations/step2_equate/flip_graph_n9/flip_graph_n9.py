"""
================================================================================
flip_graph_n9.py  --  Step-2 flip-graph connectivity check on n=9 triangulations
================================================================================

PURPOSE (BIG PICTURE)
---------------------
At n=9 the coupled-cluster matrix from
`computations/step4_laurent_block_analysis/n9/` has
  cluster rank = 90 = full,  cluster nullity = 0.
The full system retains nullity 3, all of which lives in the EXTERNAL
columns (32 triangulation cousins + 21 non-cluster step-1 survivors).

The 21 non-cluster survivors die by their own depth-1 cascades (already
verified). The 32 triangulation cousins should collapse to a single
common scalar c via Step-2's bare/special swap identities — completing
the n=9 locality theorem and closing the file on n=9.

This script verifies precisely that piece. We:

  Step 1.  Enumerate every triangulation of the 9-gon (Catalan number
           C_7 = 429 expected). Check that they group into 32 cyclic
           orbits.

  Step 2.  Build the Step-2 flip graph. For each triangulation T and
           each cyclic zone Z_{r, r+2} (r = 1..9):
             - special chord  S_r = X_{r, r+2}
             - bare chord     B_r = X_{r+1, r-1} (cyclic)
             - if T contains exactly one of {S_r, B_r}, swapping that
               chord for the other and checking the result is still a
               triangulation gives a Step-2 *swap edge*.
           Each edge corresponds to a bare/special exchange that
           Step-2 forces a coefficient-equality (up to sign) on.

  Step 3.  Lift the swap edges from triangulation-level to ORBIT-level
           (the 32 cyclic orbits) and compute connected components.

  Step 4.  Output:
             - 32 orbit reps with sizes
             - the full triangulation-level edge list
             - the orbit-level edge list (orbit-class swaps)
             - component count and component representatives

The expected outcome (per the refined Conjecture 10) is a SINGLE
connected component spanning all 32 orbits. If so, all triangulation
coefficients equal a common scalar c, and combined with the
cluster's rank-90 kill at depth-1 plus the 21 non-cluster cascades,
the n=9 locality theorem is fully proven.

OUTPUTS
-------
  outputs/n9_triangulation_orbits.json        -- 32 orbits + reps + sizes
  outputs/n9_step2_bare_swap_pairs.json       -- swap edges (tri-level + orbit-level)
  outputs/n9_flip_graph_connectivity.txt      -- human-readable verdict

CODE STYLE
----------
Per repo contribution standards (root README §7), every function has a
docstring with a LOGIC section (algorithm) and a PHYSICS / MATHEMATICS
section (what it computes in proof language).
"""

import json
import os
import sys
import time
from collections import defaultdict
from itertools import combinations

# Pull machinery from the cascade library (parametric in n).
HERE = os.path.dirname(os.path.abspath(__file__))
N8_DIR = os.path.normpath(
    os.path.join(HERE, "..", "..", "step4_laurent_block_analysis", "n8",
                 "scripts"))
sys.path.insert(0, N8_DIR)
from cascade_kill_n8 import (  # noqa: E402
    all_chords, zone_structure, is_triangulation, format_multiset,
)
from analyze_recipes import canonical_orbit_rep  # noqa: E402


# ============================================================================
# 1.  ENUMERATE TRIANGULATIONS OF THE 9-GON
# ============================================================================

def enumerate_triangulations(n=9):
    """
    Enumerate every triangulation of the n-gon as a tuple of chord pairs.

    LOGIC:  iterate every (n-3)-element subset of all_chords(n); keep the
            ones for which is_triangulation returns True.

    PHYSICS:  triangulations of the n-gon are exactly the maximal sets of
              n-3 pairwise non-crossing chords; their count is the
              Catalan number C_{n-2}.  At n=9 this gives C_7 = 429.
    """
    chords = all_chords(n)
    out = []
    for combo in combinations(chords, n - 3):
        if is_triangulation(list(combo), n):
            out.append(tuple(sorted(combo)))
    return out


def cyclic_orbit_decomposition(triangulations, n=9):
    """
    Group triangulations into Z_n cyclic orbits keyed by canonical rep.

    LOGIC:  for each triangulation, compute canonical_orbit_rep, which
            is the lex-smallest cyclic shift; group by that key.

    PHYSICS:  the cyclic Z_n action on chord sets respects all
              kill-mechanism identities, so orbit-level coefficients are
              well-defined (a_T = a_{sigma T} for any cyclic shift sigma).
              At n=9 we expect exactly 32 orbits (a known combinatorial
              count).
    """
    orbits = defaultdict(list)
    for T in triangulations:
        rep = canonical_orbit_rep(T, n)
        orbits[rep].append(T)
    return orbits


# ============================================================================
# 2.  STEP-2 BARE/SPECIAL SWAP EDGES
# ============================================================================

def step2_swap_edges(triangulations, n=9):
    """
    For each triangulation T and each zone Z_{r, r+2}, attempt the
    bare/special swap and record an edge if the swap produces another
    triangulation.

    LOGIC:
      For r in 1..n:
        Get (special S_r, bare B_r) from zone_structure(r, n).
        For each triangulation T:
          - If T contains S_r and not B_r:
              T' = (T \\ {S_r}) cup {B_r}; if T' is a triangulation,
              add edge (T, T') to the edge set (unordered pair).
          - Symmetric case: bare in T but not special.
          - If T contains BOTH or NEITHER, no Step-2 swap is defined at
            this zone for T.

    PHYSICS:  Step-2 of the kill mechanism uses the bare relation
              X_{r+1, r-1} = -X_{r, r+2} to glue pairs of multisets
              that differ only by one chord-swap.  For triangulations
              the resulting identity reads
                  a_T = (sign) * a_{T'}
              where T' = T after the swap.  Connectivity in the
              resulting "flip graph" (modulo signs) is what propagates
              equality of coefficients across all triangulations.

    Returns:
      edges_tri:   set of frozensets {T1, T2} of swap-related triangulations.
      zone_per_edge: dict edge -> list of zones at which the swap fires
                    (a single edge may be produced by multiple zones).
    """
    triangulations_set = set(triangulations)
    edges_tri = set()
    zone_per_edge = defaultdict(list)
    for r in range(1, n + 1):
        special, bare, _ = zone_structure(r, n)
        for T in triangulations:
            T_set = set(T)
            has_S = special in T_set
            has_B = bare in T_set
            if has_S and not has_B:
                T_prime_set = (T_set - {special}) | {bare}
                T_prime = tuple(sorted(T_prime_set))
                if T_prime in triangulations_set:
                    edge = frozenset([T, T_prime])
                    edges_tri.add(edge)
                    zone_per_edge[edge].append(r)
            elif has_B and not has_S:
                T_prime_set = (T_set - {bare}) | {special}
                T_prime = tuple(sorted(T_prime_set))
                if T_prime in triangulations_set:
                    edge = frozenset([T, T_prime])
                    edges_tri.add(edge)
                    zone_per_edge[edge].append(r)
    return edges_tri, dict(zone_per_edge)


def lift_edges_to_orbits(edges_tri, n=9):
    """
    Project each triangulation-level edge to an orbit-level edge.

    LOGIC:  for each edge {T1, T2}, compute canonical_orbit_rep of
            both endpoints; the orbit-level edge is the (unordered)
            pair {rep1, rep2}.  Self-loops (rep1 == rep2) are recorded
            separately.

    PHYSICS:  orbit-level connectivity is the right level to test --
              every member of an orbit shares the same coefficient by
              cyclic equivariance, so we only need orbit-level
              connectivity for the unitarity argument.
    """
    edges_orbit = set()
    self_loops = set()
    for edge in edges_tri:
        T1, T2 = list(edge)
        rep1 = canonical_orbit_rep(T1, n)
        rep2 = canonical_orbit_rep(T2, n)
        if rep1 == rep2:
            self_loops.add(rep1)
        else:
            edges_orbit.add(frozenset([rep1, rep2]))
    return edges_orbit, self_loops


# ============================================================================
# 3.  CONNECTED-COMPONENT ANALYSIS (UNION-FIND)
# ============================================================================

def connected_components(nodes, edges):
    """
    Compute connected components of a graph via union-find.

    LOGIC:  classic union-find with path compression. Each node starts
            in its own component; merge components on each edge.  Return
            list of components (each a sorted list of nodes).

    PHYSICS:  components are the maximal sets of orbits whose
              coefficients are forced equal (modulo signs) by Step-2
              alone.  A single component spanning all orbits ⇔ unitarity
              by Step-2.
    """
    parent = {node: node for node in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for edge in edges:
        a, b = list(edge)
        union(a, b)
    comps = defaultdict(list)
    for node in nodes:
        comps[find(node)].append(node)
    return list(comps.values())


# ============================================================================
# 4.  DRIVER
# ============================================================================

def main():
    n = 9
    out_dir = os.path.join(HERE, "outputs")
    os.makedirs(out_dir, exist_ok=True)

    print(f"Step-2 flip-graph connectivity check at n={n}\n")

    # Step 1: triangulations + orbits.
    print("[1/4] Enumerating triangulations...")
    t0 = time.time()
    triangulations = enumerate_triangulations(n)
    print(f"      {len(triangulations)} triangulations in "
          f"{time.time()-t0:.1f} s (expected C_7 = 429).")

    print("[2/4] Cyclic orbit decomposition...")
    orbits = cyclic_orbit_decomposition(triangulations, n)
    rep_keys_sorted = sorted(orbits.keys())
    print(f"      {len(orbits)} cyclic orbits "
          f"(expected 32 for n=9).")
    sizes = sorted({len(v) for v in orbits.values()})
    print(f"      orbit sizes: {sizes}")

    # Step 2: Step-2 swap edges.
    print("[3/4] Computing Step-2 bare/special swap edges...")
    t0 = time.time()
    edges_tri, zone_per_edge = step2_swap_edges(triangulations, n)
    edges_orbit, self_loops = lift_edges_to_orbits(edges_tri, n)
    print(f"      {len(edges_tri)} triangulation-level swap edges in "
          f"{time.time()-t0:.1f} s.")
    print(f"      {len(edges_orbit)} orbit-level edges; "
          f"{len(self_loops)} self-loops on orbits.")

    # Step 3: connected components on orbit-level graph.
    print("[4/4] Computing connected components...")
    comps = connected_components(rep_keys_sorted, edges_orbit)
    comps_sorted = sorted(comps, key=lambda c: -len(c))
    print(f"      {len(comps)} connected component(s).")
    for i, c in enumerate(comps_sorted):
        print(f"        component {i+1}: size {len(c)}")

    # ----- write outputs -----

    # Orbits manifest
    orbits_manifest = []
    for i, rep in enumerate(rep_keys_sorted, start=1):
        orbits_manifest.append({
            "orbit_id": i,
            "size": len(orbits[rep]),
            "representative": [list(c) for c in rep],
        })
    with open(os.path.join(out_dir, "n9_triangulation_orbits.json"),
              "w") as f:
        json.dump({"n": n, "n_triangulations": len(triangulations),
                   "n_orbits": len(orbits), "orbits": orbits_manifest},
                  f, indent=2)

    # Edge list
    rep_to_id = {rep: i + 1 for i, rep in enumerate(rep_keys_sorted)}
    tri_edges_serial = []
    for edge in edges_tri:
        T1, T2 = list(edge)
        rep1 = canonical_orbit_rep(T1, n)
        rep2 = canonical_orbit_rep(T2, n)
        zones = zone_per_edge[edge]
        tri_edges_serial.append({
            "T1": [list(c) for c in T1],
            "T2": [list(c) for c in T2],
            "T1_orbit": rep_to_id[rep1],
            "T2_orbit": rep_to_id[rep2],
            "zones": sorted(zones),
        })
    orbit_edges_serial = []
    for edge in edges_orbit:
        a, b = sorted(list(edge))
        orbit_edges_serial.append({
            "orbit_pair": sorted([rep_to_id[a], rep_to_id[b]]),
            "rep_a": [list(c) for c in a],
            "rep_b": [list(c) for c in b],
        })
    self_loops_serial = [
        {"orbit_id": rep_to_id[s], "rep": [list(c) for c in s]}
        for s in self_loops
    ]
    with open(os.path.join(out_dir, "n9_step2_bare_swap_pairs.json"),
              "w") as f:
        json.dump({
            "n": n,
            "n_triangulation_edges": len(edges_tri),
            "n_orbit_edges": len(edges_orbit),
            "n_self_loops": len(self_loops),
            "triangulation_edges": tri_edges_serial,
            "orbit_edges": sorted(orbit_edges_serial,
                                   key=lambda e: tuple(e["orbit_pair"])),
            "self_loops": self_loops_serial,
        }, f, indent=2)

    # Connectivity verdict.
    lines = []
    lines.append(f"Step-2 flip-graph connectivity at n={n}")
    lines.append("=" * 78)
    lines.append("")
    lines.append(f"Triangulations of the 9-gon: {len(triangulations)} "
                 f"(expected C_7 = 429)")
    lines.append(f"Cyclic orbits:               {len(orbits)} "
                 f"(expected 32)")
    lines.append(f"Orbit-size distribution:     "
                 f"{dict((sz, sum(1 for v in orbits.values() if len(v) == sz)) for sz in sorted(sizes))}")
    lines.append("")
    lines.append(f"Triangulation-level swap edges: {len(edges_tri)}")
    lines.append(f"Orbit-level swap edges:         {len(edges_orbit)}")
    lines.append(f"Orbit self-loops (swap stays in same orbit): "
                 f"{len(self_loops)}")
    lines.append("")
    lines.append(f"Connected components (orbit-level): {len(comps)}")
    if len(comps) == 1:
        lines.append("")
        lines.append("VERDICT: one connected component spanning all "
                     f"{len(orbits)} orbits.")
        lines.append("All triangulation coefficients are equated by Step-2.")
        lines.append("Combined with the cluster's rank-90 kill at depth-1")
        lines.append("and the 21 non-cluster cascades, n=9 locality is")
        lines.append("FULLY VERIFIED.")
    else:
        lines.append("")
        lines.append(f"VERDICT: {len(comps)} components -- Step-2 alone "
                     f"does NOT bridge all orbits at n=9.")
        lines.append("")
        lines.append("Component sizes (largest first):")
        for i, c in enumerate(comps_sorted):
            lines.append(f"  component {i+1}: {len(c)} orbits")
            sample_orbits = sorted(rep_to_id[rep] for rep in c)[:8]
            lines.append(f"    sample orbit IDs: {sample_orbits}")
            sample_rep = c[0]
            lines.append(f"    representative {rep_to_id[sample_rep]}: "
                         f"{format_multiset(sample_rep)}")
        lines.append("")
        lines.append("Investigate why Step-2 fails to bridge the components")
        lines.append("(may need Step-3's external columns to bridge them)")
    with open(os.path.join(out_dir, "n9_flip_graph_connectivity.txt"),
              "w") as f:
        f.write("\n".join(lines) + "\n")

    print()
    print("Outputs written to:", out_dir)
    print("  n9_triangulation_orbits.json")
    print("  n9_step2_bare_swap_pairs.json")
    print("  n9_flip_graph_connectivity.txt")


if __name__ == "__main__":
    main()
