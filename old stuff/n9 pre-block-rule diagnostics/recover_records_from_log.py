"""
Recover orbit records from a live-output log when the JSON-write was
buggy (records.append was missing).

LOGIC:  parse the log file (bm8fz93k6.output) for orbit blocks of the
        form
            Orbit N/113  (size S): M_rep = {...}
            ...
            Kill zone Z = Z_{r,r+2};  Laurent order k = K;
            Substitute U = (...);  Companion Y_U = (...);  Y_U in M_rep? = ...
            Fingerprint = ...
            Fish term scalar = ...; Cousins (all step-1 killable) = ...
            (T s)
        Or for timeouts:
            Orbit N/113  (size S): M_rep = {...}
            TIMED OUT after T s ...

        Write the recovered records into the existing JSON (merging by
        orbit_id; we don't overwrite already-present orbits).

PHYSICS:  no new computation; just data salvage.
"""
import json
import os
import re
import sys

LOG = sys.argv[1] if len(sys.argv) > 1 else None
_HERE = os.path.dirname(os.path.abspath(__file__))
JSON_PATH = os.path.normpath(os.path.join(
    _HERE, "..", "outputs", "results_cascade_n9_reps.json"))


def parse_log(log_path):
    """Yield record dicts from the log."""
    with open(log_path) as f:
        text = f.read()
    # Split on "Orbit N/" blocks.
    blocks = re.split(r"\n-{20,}\nOrbit (\d+)/(\d+)\s*", text)
    # blocks[0] = preamble; then alternating (orbit_id, total, body)
    out = []
    i = 1
    while i < len(blocks):
        try:
            orbit_id = int(blocks[i])
            i += 1
            total = int(blocks[i])
            i += 1
            body = blocks[i]
            i += 1
        except (IndexError, ValueError):
            break
        m_size = re.search(r"\(size (\d+)\)", body)
        m_rep = re.search(r"M_rep = \{([^\}]+)\}", body)
        if not m_size or not m_rep:
            continue
        size = int(m_size.group(1))
        chord_pairs = re.findall(r"\((\d+),(\d+)\)", m_rep.group(1))
        rep = [[int(a), int(b)] for a, b in chord_pairs]

        rec = {
            "orbit_id": orbit_id,
            "size": size,
            "representative": rep,
            "vertices_missed": None,
            "elapsed_s": None,
            "timed_out": False,
        }

        # Timeout?
        m_timeout = re.search(r"TIMED OUT after ([\d.]+)\s*s", body)
        if m_timeout:
            rec["timed_out"] = True
            rec["elapsed_s"] = float(m_timeout.group(1))
            rec["recipe_found"] = False
            out.append(rec)
            continue

        # Successful recipe?
        m_zone = re.search(
            r"Kill zone Z = Z_\{(\d+),(\d+)\};\s+Laurent order k = (\d+);",
            body)
        m_sub = re.search(
            r"Substitute U = \((\d+), (\d+)\);\s+Companion Y_U = \((\d+), (\d+)\);"
            r"\s+Y_U in M_rep\? = (True|False)",
            body)
        m_fp = re.search(r"Fingerprint = (.+)", body)
        m_scalar = re.search(
            r"Fish term scalar = (-?\d+);\s+Cousins \(all step-1 killable\) = (\d+)",
            body)
        m_t = re.search(r"\(([\d.]+) s\)", body)
        if m_zone and m_sub and m_fp and m_scalar:
            rec["recipe_found"] = True
            rec["zone_r"] = int(m_zone.group(1))
            rec["order_k"] = int(m_zone.group(3))
            rec["substitute_U"] = [int(m_sub.group(1)), int(m_sub.group(2))]
            rec["companion_Y_U"] = [int(m_sub.group(3)), int(m_sub.group(4))]
            rec["Y_U_in_M_rep"] = (m_sub.group(5) == "True")
            rec["fingerprint"] = m_fp.group(1).strip()
            rec["fish_scalar"] = int(m_scalar.group(1))
            rec["n_cousins"] = int(m_scalar.group(2))
            if m_t:
                rec["elapsed_s"] = float(m_t.group(1))
            out.append(rec)
            continue

        # Otherwise: incomplete / interrupted -- skip.

    return out


def main():
    if not LOG:
        print("Usage: recover_records_from_log.py <log_path>")
        sys.exit(1)
    new_records = parse_log(LOG)
    print(f"Parsed {len(new_records)} recoverable records from {LOG}.")
    for r in new_records[:5]:
        print(f"  orbit {r['orbit_id']}  found={r.get('recipe_found')}  "
              f"timed_out={r.get('timed_out')}  elapsed={r.get('elapsed_s')}")

    with open(JSON_PATH) as f:
        existing = json.load(f)
    existing_records = existing.get("records", [])
    existing_ids = {r["orbit_id"] for r in existing_records}
    merged = list(existing_records)
    added = 0
    for r in new_records:
        if r["orbit_id"] not in existing_ids:
            merged.append(r)
            added += 1
    print(f"Adding {added} new records (skipped {len(new_records) - added} "
          f"already present).")

    successes = sum(1 for r in merged if r.get("recipe_found"))
    failures_list = [r["representative"] for r in merged
                     if not r.get("recipe_found")]
    out_data = {
        "n": existing.get("n", 9),
        "orbits_total": existing.get("orbits_total", 113),
        "completed": len(merged),
        "successes": successes,
        "failures": failures_list,
        "records": merged,
    }
    with open(JSON_PATH, "w") as f:
        json.dump(out_data, f, indent=2)
    print(f"Wrote {len(merged)} records to {JSON_PATH}.")


if __name__ == "__main__":
    main()
