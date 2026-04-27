# Diagnostic for orbit 28

M_rep = {(1,3),(1,4),(1,8),(2,4),(4,6),(6,8)}

V_used = [1, 2, 3, 4, 6, 8]
V_missed = [5, 7, 9]

## Step-1 (K1)+(K2) status across all zones

| zone | special | bare | ell_Z | K1 fail (special in M) | K1 fail (bare in M) | K2 fails (pairs) | killable? |
|---|---|---|---|---|---|---|---|
| Z_{1,3} | (1, 3) | (2, 9) | 2 | YES | no | ((1, 4),(2, 4)) | NOT killable |
| Z_{2,4} | (2, 4) | (1, 3) | 2 | YES | YES | — | NOT killable |
| Z_{3,5} | (3, 5) | (2, 4) | 3 | no | YES | ((1, 3),(1, 4)) | NOT killable |
| Z_{4,6} | (4, 6) | (3, 5) | 1 | YES | no | — | NOT killable |
| Z_{5,7} | (5, 7) | (4, 6) | 2 | no | YES | — | NOT killable |
| Z_{6,8} | (6, 8) | (5, 7) | 1 | YES | no | — | NOT killable |
| Z_{7,9} | (7, 9) | (6, 8) | 2 | no | YES | — | NOT killable |
| Z_{8,1} | (1, 8) | (7, 9) | 1 | YES | no | — | NOT killable |
| Z_{9,2} | (2, 9) | (1, 8) | 3 | no | YES | — | NOT killable |

Confirmed: M_rep is a genuine step-1 survivor (not killable at any of the 9 zones).

## Depth-1 cascade attempts

Tried depth-1 cascades across all zones in 318.4 s; found 6 fingerprint equations containing M_rep with nonzero scalar.
Of those, 2 have ALL cousins step-1 killable (= valid kill recipe).

### Valid recipe @ depth 1, zone Z_{5,7}, order 3

- Fingerprint: `X_5_8/(X_1_3*X_1_4*X_1_8*X_2_4)`
- Chosen substitutes (with multiplicity): [([6, 8], 1)]
- M_rep scalar: 1
- Cousins (7): all step-1 killable.
  - [[1, 3], [1, 4], [1, 6], [1, 8], [2, 4], [6, 8]]: scalar=1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [2, 6], [6, 8]]: scalar=1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [3, 6], [6, 8]]: scalar=1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [5, 7], [6, 8]]: scalar=-1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [5, 8], [6, 8]]: scalar=-1, kill at Z_{4,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [6, 8], [6, 8]]: scalar=2, kill at Z_{4,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [6, 8], [6, 9]]: scalar=1, kill at Z_{4,*}

### Valid recipe @ depth 1, zone Z_{7,9}, order 3

- Fingerprint: `X_1_7/(X_1_3*X_1_4*X_2_4*X_4_6)`
- Chosen substitutes (with multiplicity): [([1, 8], 1)]
- M_rep scalar: 1
- Cousins (7): all step-1 killable.
  - [[1, 3], [1, 4], [1, 7], [1, 8], [2, 4], [4, 6]]: scalar=-1, kill at Z_{6,*}
  - [[1, 3], [1, 4], [1, 8], [1, 8], [2, 4], [4, 6]]: scalar=2, kill at Z_{6,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [2, 8], [4, 6]]: scalar=1, kill at Z_{6,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [3, 8], [4, 6]]: scalar=1, kill at Z_{6,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [4, 6], [4, 8]]: scalar=1, kill at Z_{6,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [4, 6], [5, 8]]: scalar=1, kill at Z_{6,*}
  - [[1, 3], [1, 4], [1, 8], [2, 4], [4, 6], [7, 9]]: scalar=-1, kill at Z_{6,*}

=== STOP: depth-1 closes orbit 28. ===

