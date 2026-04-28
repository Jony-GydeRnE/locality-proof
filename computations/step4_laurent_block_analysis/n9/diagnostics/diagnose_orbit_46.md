# Diagnostic for orbit 46

M_rep = {(1,3),(1,4),(1,8),(4,6),(5,7),(7,9)}

V_used = [1, 3, 4, 5, 6, 7, 8, 9]
V_missed = [2]

## Step-1 (K1)+(K2) status across all zones

| zone | special | bare | ell_Z | K1 fail (special in M) | K1 fail (bare in M) | K2 fails (pairs) | killable? |
|---|---|---|---|---|---|---|---|
| Z_{1,3} | (1, 3) | (2, 9) | 1 | YES | no | — | NOT killable |
| Z_{2,4} | (2, 4) | (1, 3) | 1 | no | YES | — | NOT killable |
| Z_{3,5} | (3, 5) | (2, 4) | 2 | no | no | ((1, 3),(1, 4)) | NOT killable |
| Z_{4,6} | (4, 6) | (3, 5) | 2 | YES | no | — | NOT killable |
| Z_{5,7} | (5, 7) | (4, 6) | 2 | YES | YES | — | NOT killable |
| Z_{6,8} | (6, 8) | (5, 7) | 2 | no | YES | — | NOT killable |
| Z_{7,9} | (7, 9) | (6, 8) | 2 | YES | no | — | NOT killable |
| Z_{8,1} | (1, 8) | (7, 9) | 2 | YES | YES | — | NOT killable |
| Z_{9,2} | (2, 9) | (1, 8) | 3 | no | YES | — | NOT killable |

Confirmed: M_rep is a genuine step-1 survivor (not killable at any of the 9 zones).

## Depth-1 cascade attempts

Tried depth-1 cascades across all zones in 236.1 s; found 7 fingerprint equations containing M_rep with nonzero scalar.
Of those, 0 have ALL cousins step-1 killable (= valid kill recipe).

Top 5 partial fingerprint equations at depth 1 (cousins not all step-1 killable):

- Z_{3,5}, fp `X_3_6/(X_1_3*X_1_8*X_5_7*X_7_9)`: M scalar=1, 3 survivor-cousins out of 7
  - survivor cousin: [[1, 3], [1, 8], [2, 4], [4, 6], [5, 7], [7, 9]] (scalar 1)
  - survivor cousin: [[1, 3], [1, 8], [3, 5], [4, 6], [5, 7], [7, 9]] (scalar -1)
  - survivor cousin: [[1, 3], [1, 8], [3, 6], [4, 6], [5, 7], [7, 9]] (scalar -1)

- Z_{3,5}, fp `1/(X_1_8*X_5_7*X_7_9)`: M scalar=1, 5 survivor-cousins out of 123
  - survivor cousin: [[1, 3], [1, 4], [1, 8], [3, 5], [5, 7], [7, 9]] (scalar -1)
  - survivor cousin: [[1, 3], [1, 4], [1, 8], [4, 7], [5, 7], [7, 9]] (scalar 1)
  - survivor cousin: [[1, 4], [1, 8], [2, 4], [3, 5], [5, 7], [7, 9]] (scalar 1)

- Z_{4,6}, fp `X_4_7/(X_1_3*X_1_4*X_1_8*X_7_9)`: M scalar=-1, 3 survivor-cousins out of 7
  - survivor cousin: [[1, 3], [1, 4], [1, 5], [1, 8], [5, 7], [7, 9]] (scalar 1)
  - survivor cousin: [[1, 3], [1, 4], [1, 8], [3, 5], [5, 7], [7, 9]] (scalar 1)
  - survivor cousin: [[1, 3], [1, 4], [1, 8], [4, 7], [5, 7], [7, 9]] (scalar -1)

- Z_{6,8}, fp `X_6_9/(X_1_3*X_1_4*X_1_8*X_4_6)`: M scalar=1, 3 survivor-cousins out of 7
  - survivor cousin: [[1, 3], [1, 4], [1, 8], [4, 6], [4, 7], [7, 9]] (scalar 1)
  - survivor cousin: [[1, 3], [1, 4], [1, 8], [4, 6], [6, 8], [7, 9]] (scalar -1)
  - survivor cousin: [[1, 3], [1, 4], [1, 8], [4, 6], [6, 9], [7, 9]] (scalar -1)

- Z_{7,9}, fp `X_1_7/(X_1_3*X_1_4*X_4_6*X_5_7)`: M scalar=-1, 3 survivor-cousins out of 7
  - survivor cousin: [[1, 3], [1, 4], [1, 7], [1, 8], [4, 6], [5, 7]] (scalar -1)
  - survivor cousin: [[1, 3], [1, 4], [1, 8], [4, 6], [5, 7], [5, 8]] (scalar 1)
  - survivor cousin: [[1, 3], [1, 4], [1, 8], [4, 6], [5, 7], [6, 8]] (scalar 1)

## Depth-2 cascade attempts

Tried depth-2 cascades across all zones in 175.4 s; found 2 fingerprint equations containing M_rep with nonzero scalar.
Of those, 0 have ALL cousins step-1 killable (= valid kill recipe).

Top 2 partial fingerprint equations at depth 2 (cousins not all step-1 killable):

- Z_{3,5}, fp `X_3_6/(X_1_8*X_5_7*X_7_9)`: M scalar=1, 1 survivor-cousins out of 39
  - survivor cousin: [[1, 4], [1, 8], [2, 4], [4, 6], [5, 7], [7, 9]] (scalar -1)

- Z_{9,2}, fp `X_3_9*X_4_9/(X_4_6*X_5_7*X_7_9)`: M scalar=-1, 4 survivor-cousins out of 8
  - survivor cousin: [[1, 3], [1, 4], [1, 7], [4, 6], [5, 7], [7, 9]] (scalar -1)
  - survivor cousin: [[1, 3], [1, 4], [2, 9], [4, 6], [5, 7], [7, 9]] (scalar 1)
  - survivor cousin: [[1, 3], [1, 4], [3, 9], [4, 6], [5, 7], [7, 9]] (scalar 1)

## Depth-3 cascade attempts

Tried depth-3 cascades across all zones in 0.0 s; found 0 fingerprint equations containing M_rep with nonzero scalar.
Of those, 0 have ALL cousins step-1 killable (= valid kill recipe).

## Depth-4 cascade attempts

Tried depth-4 cascades across all zones in 0.0 s; found 0 fingerprint equations containing M_rep with nonzero scalar.
Of those, 0 have ALL cousins step-1 killable (= valid kill recipe).


=== UNRESOLVED at depth <= 4. ===

