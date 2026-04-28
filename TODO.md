# TODO — pickup notes for future sessions

Last touched: 2026-04-28.

## Current state of the repo

- **n=9 locality is fully proven** from the cyclic 1-zero conditions alone.
  Consolidated artifact + verifier:
  `computations/step4_laurent_block_analysis/n9/outputs/n9_locality_status.md`
  (verdict: 113/113 orbit reps accounted for, 23 single-orbit cascades
  + 90-orbit block-rule, all ✓).
- **n=9 unitarity** (all triangulation $c$-values equal one common $c$)
  is handled by the d-subset argument in
  `paper/Details missing in Hidden zero -> unitarity Rodina proof/details for Rodina D subset argument/`
  and notes 11, 18 — independent of the cascade machinery.
- The 4 "untouched" Step-2 components at $n=9$ are TRIANGULATION
  components (chord-length signature `(2,2,3,3,4,4)`, fan-class). They
  are LOCAL terms, not a locality gap. Equating them is the
  unitarity question above.
- The latest paper write-up is
  `paper/step3 block rule/step3_block_rule_n9_REVISED.{tex,pdf}` and the
  matching README `step3_block_rule_README_REVISED.md`.

## Pending — big repo reorganization (Notion `📋 Claude Code Prompt`)

**Skipped this session intentionally.** A future session can execute
this if/when the user wants. The Notion page
`📋 Claude Code Prompt: Locality-Proof Repo Reorganization (5 Changes)`
under the master plan dashboard has the full spec. Summary:

1. **Merge `paper/` and `computations/` into unified `stepN-…/` folders**,
   each with `paper/`, `code/`, `outputs/` subfolders. Top-level layout:
   ```
   step0-geometric-foundations/
   step1-layer0-kill/
   step2-equate-unitarity/
   step3-laurent-cascade/
   step4-block-rule/
   foundations/                 ← was paper/Details missing…
   full-nullspace-verification/
   notes/, existing-literature/, old-stuff/   ← rename hyphenated
   ```
   Use `git mv`. Delete empty `computations/` and `paper/` after.
2. Update root README §3 status table and add an executive summary
   after the §6 master plan table.
3. Replace "core paper does not yet exist" wording in `paper/README.md`.
4. Add `notes/README.md` with explicit "load-bearing notes 12-18"
   guidance + a brief reference / dated section.
5. Re-execute the moves end-to-end with the verification checklist
   from the Notion prompt (every .py reachable, every .pdf in
   `paper/`, READMEs cross-linked, no broken paths, `git status`
   shows only renames).

**Risks for the future session:**

- Many n9 scripts import from `../../n8/scripts/cascade_kill_n8.py`.
  After the move those imports will break and need rewiring (it's the
  same problem from the previous reorg; pattern: `OUT_DIR = ../outputs`,
  `N8_DIR = ../../n8/scripts` → updated to the new layout).
- READMEs at root, `paper/`, every stepN/, and several scripts
  hardcode paths like `../outputs/orbits_n9.json`. All need
  updating.
- Folders with spaces (`paper/step1 Kill technique and statistics/`,
  `paper/Details missing in Hidden zero -> unitarity Rodina proof/`)
  need careful quoting in `git mv`.

## Smaller follow-ups, not blocking

- The `paper/step3 block rule/step3_block_rule_n9_REVISED.{tex,pdf,md}`
  files supersede the unREVISED versions; once the user is happy
  with REVISED, delete the unREVISED originals (or rename REVISED →
  canonical).
- `notes/18. d subset rigorous proof.pdf` is now load-bearing for
  unitarity — make sure the d-subset paper in
  `paper/Details missing.../details for Rodina D subset argument/`
  cites it explicitly.
- Extend the cluster-matrix verification to $n=10$ (same
  block-rule kill should work; not yet attempted).

## What NOT to do

- Don't touch `notes/`, `existing literature/`, `old stuff/` contents —
  only rename folders to hyphenated form if executing CHANGE 1.
- Don't run another full $n=9$ cascade or cluster matrix — those are
  done and recorded; pull from `outputs/`.
- Don't rebuild the Step-2 flip graph — `flip_graph_n9` is final.
- Don't try to "kill" the 4 fan-class triangulation components —
  those are LOCAL and survive by design.
