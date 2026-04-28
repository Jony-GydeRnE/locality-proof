# Diagnostic for orbit 68

M_rep = {(1,3),(1,4),(4,6),(4,9),(6,8),(7,9)}

V_used = [1, 3, 4, 6, 7, 8, 9]
V_missed = [2, 5]

## Step-1 (K1)+(K2) status across all zones

| zone | special | bare | ell_Z | K1 fail (special in M) | K1 fail (bare in M) | K2 fails (pairs) | killable? |
|---|---|---|---|---|---|---|---|
| Z_{1,3} | (1, 3) | (2, 9) | 1 | YES | no | — | NOT killable |
| Z_{2,4} | (2, 4) | (1, 3) | 1 | no | YES | — | NOT killable |
| Z_{3,5} | (3, 5) | (2, 4) | 3 | no | no | ((1, 3),(1, 4)) | NOT killable |
| Z_{4,6} | (4, 6) | (3, 5) | 1 | YES | no | — | NOT killable |
| Z_{5,7} | (5, 7) | (4, 6) | 2 | no | YES | — | NOT killable |
| Z_{6,8} | (6, 8) | (5, 7) | 2 | YES | no | — | NOT killable |
| Z_{7,9} | (7, 9) | (6, 8) | 2 | YES | YES | — | NOT killable |
| Z_{8,1} | (1, 8) | (7, 9) | 2 | no | YES | — | NOT killable |
| Z_{9,2} | (2, 9) | (1, 8) | 2 | no | no | ((4, 9),(1, 4)) | NOT killable |

Confirmed: M_rep is a genuine step-1 survivor (not killable at any of the 9 zones).

## Depth-1 cascade attempts

Tried depth-1 cascades across all zones in 224.4 s; found 7 fingerprint equations containing M_rep with nonzero scalar.
Of those, 1 have ALL cousins step-1 killable (= valid kill recipe).

### Valid recipe @ depth 1, zone Z_{5,7}, order 3

- Fingerprint: `X_5_8/(X_1_3*X_1_4*X_4_9*X_7_9)`
- Chosen substitutes (with multiplicity): [([6, 8], 1)]
- M_rep scalar: 1
- Cousins (7): all step-1 killable.
  - [[1, 3], [1, 4], [1, 6], [4, 9], [6, 8], [7, 9]]: scalar=1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [2, 6], [4, 9], [6, 8], [7, 9]]: scalar=1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [3, 6], [4, 9], [6, 8], [7, 9]]: scalar=1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [4, 9], [5, 7], [6, 8], [7, 9]]: scalar=-1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [4, 9], [5, 8], [6, 8], [7, 9]]: scalar=-1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [4, 9], [6, 8], [6, 8], [7, 9]]: scalar=2, kill at Z_{4,*}
  - [[1, 3], [1, 4], [4, 9], [6, 8], [6, 9], [7, 9]]: scalar=1, kill at Z_{4,*}

=== STOP: depth-1 closes orbit 68. ===

