"""
================================================================================
cascade_kill_orbit_reps_n9.py  --  Run the depth-1 Laurent cascade only on
                                     the 113 cyclic-orbit reps at n=9
================================================================================

PURPOSE (BIG PICTURE)
---------------------
At n=9 there are 1 011 non-triangulation step-1 survivors but only 113
cyclic-orbit equivalence classes (one of size 3, the others all size 9).
The cascade kill mechanism is cyclically equivariant -- a recipe found
for one orbit member determines the recipes for all members of that
orbit by cyclic shift.  So instead of running the cascade 1 011 times,
we only need to run it 113 times: once per orbit representative.

This script:

  (1) Loads the orbit manifest from `orbits_n9.json`
      (produced by `orbit_decomposition.py`).
  (2) For each representative, runs `search_depth1_recipe` (imported
      from `../cascade_n8/cascade_kill_n8.py` -- it is parametric in n).
  (3) For each successful representative, logs:
          orbit size, M_rep, kill recipe (Z, k, fingerprint, companion U,
          Y_U), whether Y_U appears in M_rep at the chosen Z.

OUTPUTS
-------
  results_cascade_n9_reps.txt  -- human-readable trace, one block per rep.
  results_cascade_n9_reps.json -- machine-readable record per rep.

CODE STYLE
----------
Per repo contribution standards: each function has a LOGIC docstring
section plus a PHYSICS / MATHEMATICS section.
"""

import json
import multiprocessing
import os
import sys
import time
import sympy as sp

# Pull the cascade machinery from cascade_n8/ (parametric in n).
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.normpath(os.path.join(HERE, "..", "outputs"))
DIAG_DIR = os.path.normpath(os.path.join(HERE, "..", "diagnostics"))
N8_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "n8", "scripts"))
sys.path.insert(0, N8_DIR)

from cascade_kill_n8 import (  # noqa: E402
    zone_structure, search_depth1_recipe, format_multiset,
    fish_substitutes_at_zone,
)


def companion_for_substitute(substitute, zone_r, n):
    """
    Look up the companion of `substitute` on zone Z_{zone_r, ...}.

    LOGIC:  scan the (companion, substitute) pairs of zone_structure
            for the matching substitute; return its companion.

    PHYSICS:  the companion Y_U of substitute U at zone Z is the
              row-r chord with the same far endpoint as U; the depth-1
              fingerprint X_{Y_U} / (rest of M's free factors) is
              produced by the Laurent tail of 1/(Y_U - X_S).
    """
    _, _, pairs = zone_structure(zone_r, n)
    for comp, sub in pairs:
        if sub == substitute:
            return comp
    return None


def _cascade_worker(orbit_record, n, max_extra_orders, queue):
    """
    Subprocess worker: run search_depth1_recipe and put the result on
    the queue.

    LOGIC:  imports here so the import side-effects (sympy startup) are
            paid in the worker, not the parent.

    PHYSICS:  isolates each rep's cascade in its own process so we can
              terminate it cleanly if it exceeds the per-rep timeout.
    """
    try:
        # Re-import inside the worker (fresh sympy state).
        import sys
        sys.path.insert(0, N8_DIR)
        from cascade_kill_n8 import search_depth1_recipe
        M_rep_list = orbit_record["representative"]
        M_rep = tuple(tuple(c) for c in M_rep_list)
        rec = search_depth1_recipe(M_rep, n=n,
                                   max_extra_orders=max_extra_orders,
                                   verbose=False)
        # Convert to a serialisable form.
        if rec is None:
            queue.put(None)
        else:
            ser = {
                "zone_r": rec["zone_r"],
                "order": rec["order"],
                "fingerprint": str(sp.simplify(rec["fingerprint"])),
                "companion": list(rec["companion"]),
                "fish_scalar": int(rec["fish_scalar"]),
                "n_cousins": len(rec["cousin_kills"]),
            }
            queue.put(ser)
    except Exception as e:
        queue.put({"_error": str(e)})


def run_one_rep(orbit_record, n=9, max_extra_orders=2,
                timeout_seconds=600):
    """
    Run the depth-1 cascade search on one orbit representative.

    LOGIC:  call search_depth1_recipe (from cascade_n8/cascade_kill_n8.py);
            if found, identify the substitute U used and look up its
            companion Y_U; check if Y_U is in M_rep at the chosen zone.

    PHYSICS:  one search per orbit rep; cyclic equivariance lifts the
              result to the whole orbit.
    """
    """
    Run the depth-1 cascade search on one orbit representative, with a
    per-rep timeout enforced via a subprocess.

    LOGIC:  spawn a worker process; wait up to `timeout_seconds`; if
            the worker exceeds the budget, terminate it and mark the
            rep as `timed_out`.

    PHYSICS:  some representatives at n=9 (especially those with
              repeated chords producing dense candidate pools) take
              an order of magnitude longer than typical; the timeout
              keeps the whole-run wall-clock bounded while flagging
              those reps for separate, deeper investigation.
    """
    M_rep_list = orbit_record["representative"]
    M_rep = tuple(tuple(c) for c in M_rep_list)
    t0 = time.time()

    queue = multiprocessing.Queue()
    proc = multiprocessing.Process(
        target=_cascade_worker,
        args=(orbit_record, n, max_extra_orders, queue),
    )
    proc.start()
    proc.join(timeout=timeout_seconds)
    elapsed = time.time() - t0
    timed_out = False
    if proc.is_alive():
        proc.terminate()
        proc.join(5)
        if proc.is_alive():
            proc.kill()
        timed_out = True
        worker_result = None
    else:
        try:
            worker_result = queue.get_nowait()
        except Exception:
            worker_result = None

    out = {
        "orbit_id": orbit_record["orbit_id"],
        "size": orbit_record["size"],
        "representative": [list(c) for c in M_rep],
        "vertices_missed": orbit_record["vertices_missed"],
        "elapsed_s": elapsed,
        "timed_out": timed_out,
    }
    if timed_out or worker_result is None or worker_result.get("_error"):
        out["recipe_found"] = False
        if worker_result and worker_result.get("_error"):
            out["error"] = worker_result["_error"]
        return out

    zone_r = worker_result["zone_r"]
    companion = tuple(worker_result["companion"])
    # Identify the substitute U whose companion is `companion` at this zone.
    sub_in_rep = None
    for comp, sub in zone_structure(zone_r, n)[2]:
        if comp == companion and sub in M_rep:
            sub_in_rep = sub
            break
    Y_in_M = (companion in M_rep)

    out.update({
        "recipe_found": True,
        "zone_r": zone_r,
        "order_k": worker_result["order"],
        "fingerprint": worker_result["fingerprint"],
        "companion_Y_U": list(companion),
        "substitute_U": list(sub_in_rep) if sub_in_rep else None,
        "Y_U_in_M_rep": Y_in_M,
        "fish_scalar": worker_result["fish_scalar"],
        "n_cousins": worker_result["n_cousins"],
    })
    return out


def main(orbits_json=None, output_txt=None, output_json=None,
         limit=None, timeout_seconds=600, skip_done=True):
    """
    Top-level driver.

    LOGIC:  load orbit manifest -> iterate over reps -> run cascade per rep
            -> write trace + JSON.

    PHYSICS:  the trace is the per-orbit kill record at n=9; combined
              with the orbit manifest, it certifies (modulo the depth-1
              cascade hypothesis itself) that all 1 011 step-1 survivors
              die.
    """
    n = 9
    if orbits_json is None:
        orbits_json = os.path.join(OUT_DIR, "orbits_n9.json")
    if output_txt is None:
        output_txt = os.path.join(OUT_DIR, "results_cascade_n9_reps.txt")
    if output_json is None:
        output_json = os.path.join(OUT_DIR, "results_cascade_n9_reps.json")

    with open(orbits_json) as f:
        manifest = json.load(f)
    reps = manifest["orbits"]

    # If a previous results file exists and skip_done is True, preload its
    # records and only run reps not yet present (so resume after a kill).
    existing_records = []
    done_orbit_ids = set()
    if skip_done and os.path.exists(output_json):
        try:
            with open(output_json) as f:
                prev = json.load(f)
            existing_records = prev.get("records", [])
            done_orbit_ids = {r["orbit_id"] for r in existing_records}
            print(f"Resume: found {len(done_orbit_ids)} already-completed "
                  f"reps in {output_json}; will skip them.")
        except Exception as e:
            print(f"Could not parse existing {output_json}: {e}")

    if done_orbit_ids:
        reps = [r for r in reps if r["orbit_id"] not in done_orbit_ids]

    if limit is not None:
        reps = reps[:limit]
        print(f"NOTE: limit set; only running next {len(reps)} reps.")

    print(f"Running depth-1 cascade on {len(reps)} orbit representatives "
          f"at n={n}...")

    log_lines = []
    def log(s):
        log_lines.append(s)
        print(s)
        sys.stdout.flush()

    log("=" * 78)
    log(f"Cascade verification on n=9 cyclic-orbit representatives")
    log("=" * 78)
    log("")
    log(f"By cyclic equivariance, success on a representative implies")
    log(f"success on every member of its orbit.  At n=9 there are 1 011")
    log(f"step-1 survivors falling into 113 orbits (112 size-9, 1 size-3),")
    log(f"so 113 cascade runs cover all 1 011 cases.")
    log("")

    records = list(existing_records)
    successes = sum(1 for r in records if r.get("recipe_found"))
    failures = [r["representative"] for r in records
                if not r.get("recipe_found")]
    t_start = time.time()
    for r in reps:
        log("-" * 78)
        log(f"Orbit {r['orbit_id']}/{len(manifest['orbits'])}  "
            f"(size {r['size']}): "
            f"M_rep = {format_multiset(tuple(tuple(c) for c in r['representative']))}")
        log("-" * 78)
        rec = run_one_rep(r, n=n, max_extra_orders=2,
                          timeout_seconds=timeout_seconds)
        records.append(rec)
        if rec.get("timed_out"):
            log(f"  TIMED OUT after {rec['elapsed_s']:.1f} s "
                f"(budget {timeout_seconds} s); marked for follow-up.")
            failures.append(r["representative"])
            log("")
            with open(output_json, "w") as f:
                json.dump({"n": n, "orbits_total": len(manifest["orbits"]),
                           "completed": len(records),
                           "successes": successes,
                           "failures": failures,
                           "records": records}, f, indent=2)
            continue
        if rec["recipe_found"]:
            sp_z, _, _ = zone_structure(rec["zone_r"], n)
            log(f"  Kill zone Z = Z_{{{sp_z[0]},{sp_z[1]}}};  Laurent order "
                f"k = {rec['order_k']};")
            log(f"  Substitute U = {tuple(rec['substitute_U']) if rec['substitute_U'] else None};  "
                f"Companion Y_U = {tuple(rec['companion_Y_U'])};  "
                f"Y_U in M_rep? = {rec['Y_U_in_M_rep']}")
            log(f"  Fingerprint = {rec['fingerprint']}")
            log(f"  Fish term scalar = {rec['fish_scalar']};  "
                f"Cousins (all step-1 killable) = {rec['n_cousins']}")
            log(f"  ({rec['elapsed_s']:.1f} s)")
            successes += 1
        else:
            log(f"  !!! No depth-1 recipe found  ({rec['elapsed_s']:.1f} s).")
            failures.append(r["representative"])
        log("")
        # Save partial progress as we go (in case of interruption).
        with open(output_json, "w") as f:
            json.dump({"n": n, "orbits_total": len(reps),
                       "completed": len(records),
                       "successes": successes,
                       "failures": failures,
                       "records": records}, f, indent=2)

    total = time.time() - t_start
    log("=" * 78)
    log(f"Summary: {successes} / {len(reps)} orbit reps killed by depth-1 "
        f"cascade in {total:.1f} s ({total/60.0:.1f} min).")
    if failures:
        log(f"Failures ({len(failures)}):")
        for f in failures:
            log(f"  {f}")
    else:
        log(f"All {len(reps)} orbit reps die via depth-1 Laurent cascade.")
        total_survivors = sum(r["size"] for r in reps)
        log(f"By cyclic equivariance, this lifts to all {total_survivors}")
        log(f"step-1 survivors at n=9.")
    log("=" * 78)

    with open(output_txt, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\nResults written to {output_txt} and {output_json}.")


if __name__ == "__main__":
    # CLI:
    #   python3 cascade_kill_orbit_reps_n9.py
    #   python3 cascade_kill_orbit_reps_n9.py LIMIT
    #   python3 cascade_kill_orbit_reps_n9.py LIMIT TIMEOUT_SECONDS
    limit = None
    timeout_seconds = 300  # 5-min budget per rep by default for the resume run.
    if len(sys.argv) > 1:
        limit = int(sys.argv[1])
    if len(sys.argv) > 2:
        timeout_seconds = int(sys.argv[2])
    print(f"Per-rep timeout: {timeout_seconds} s.")
    main(limit=limit, timeout_seconds=timeout_seconds, skip_done=True)
