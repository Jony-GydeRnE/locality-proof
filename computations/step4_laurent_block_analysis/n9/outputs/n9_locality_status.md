# n=9 locality status — consolidated artifact

This file accounts line-by-line for every cyclic-orbit representative of n=9 Step-1 survivors, showing which kill mechanism eliminates its coefficient. It is the one-stop reference for the n=9 locality theorem.

## 1. Header counts

| Quantity | Value |
|---|---:|
| n | 9 |
| Total non-triangulation size-6 multisets | 905 763 |
| Step-1-killable count | 904 752 |
| Step-1 non-tri survivor count (= # non-tri multisets uncaught at every zone) | 1 011 |
| Cyclic-orbit count (Z_9 action) | 113 |
| Cluster orbits (BFS from orbit 22) | 90 |
| Non-cluster orbits | 23 |
| Cluster matrix dimensions | 506 rows × 90 cluster cols + 53 external cols |
| Cluster rank | **90** |
| Cluster nullity (= residual r) | **0** |

## 2. Per-orbit account

All 113 cyclic-orbit representatives below.  Cluster orbits are killed by the *block-rule* (cluster matrix has rank = full = 90).  Non-cluster orbits are each killed by their own single-orbit depth-1 cascade.

Locality status (`Loc`) of each rep is independently verified: it should be `CROSSING` or `DOUBLE_POLE` (both non-local) for every one of the 113 orbits, since the orbit manifest is filtered to non-tri step-1 survivors. A `TRIANGULATION` here would indicate a manifest bug.

| Orbit | Cluster? | Loc | Representative | Kill mechanism / recipe |
|---:|:---:|:---:|---|---|
| 1 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,4),(1,8),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 2 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,7),(1,8),(3,5),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 3 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,7),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 4 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,8),(2,4),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 5 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,8),(3,5),(3,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 6 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,8),(3,5),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 7 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,8),(3,5),(5,7),(5,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 8 | NON_CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,8),(3,5),(5,7),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{4,6}, k=3, U=(5, 7), Y_U=(4, 7); fp = X_4_7/(X_1_3**2*X_1_8*X_6_8); cousins = 7 (all step-1 killable) |
| 9 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,8),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 10 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,8),(3,5),(5,8),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 11 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(1,8),(3,6),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 12 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(2,9),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 13 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,3),(3,5),(3,9),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 14 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,4),(1,4),(1,8),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 15 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,5),(1,8),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{7,9}, k=3, U=(1, 8), Y_U=(1, 7); fp = X_1_7/(X_1_3*X_1_4*X_1_5*X_4_6); cousins = 7 (all step-1 killable) |
| 16 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,5),(1,8),(5,7),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 17 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,5),(1,8),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 18 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,5),(2,9),(5,7),(7,9)}` | SINGLE-ORBIT CASCADE: Z_{6,8}, k=3, U=(7, 9), Y_U=(6, 9); fp = X_6_9/(X_1_3*X_1_4*X_1_5*X_2_9); cousins = 7 (all step-1 killable) |
| 19 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,5),(3,9),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 20 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,5),(4,9),(5,7),(7,9)}` | SINGLE-ORBIT CASCADE: Z_{3,5}, k=3, U=(4, 9), Y_U=(3, 9); fp = X_3_9/(X_1_3*X_1_5*X_5_7*X_7_9); cousins = 7 (all step-1 killable) |
| 21 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,7),(1,8),(3,5),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 22 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,7),(1,8),(4,6),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 23 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,7),(1,8),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 24 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,7),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 25 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,7),(4,6),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 26 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,7),(4,6),(6,8),(7,9)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=3, U=(6, 8), Y_U=(5, 8); fp = X_5_8/(X_1_3*X_1_4*X_1_7*X_7_9); cousins = 7 (all step-1 killable) |
| 27 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,7),(4,6),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 28 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(2,4),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 29 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(2,5),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{4,6}, k=3, U=(2, 5), Y_U=(2, 4); fp = X_2_4/(X_1_3*X_1_4*X_1_8*X_6_8); cousins = 7 (all step-1 killable) |
| 30 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(2,6),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=4, U=(2, 6), Y_U=(2, 5); fp = X_2_5/(X_1_3*X_1_4*X_1_8); cousins = 39 (all step-1 killable) |
| 31 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(2,7),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=3, U=(6, 8), Y_U=(5, 8); fp = X_5_8/(X_1_3*X_1_4*X_1_8*X_2_7); cousins = 7 (all step-1 killable) |
| 32 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(2,8),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=3, U=(6, 8), Y_U=(5, 8); fp = X_5_8/(X_1_3*X_1_4*X_1_8*X_2_8); cousins = 7 (all step-1 killable) |
| 33 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(2,9),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=3, U=(6, 8), Y_U=(5, 8); fp = X_5_8/(X_1_3*X_1_4*X_1_8*X_2_9); cousins = 7 (all step-1 killable) |
| 34 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(3,5),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 35 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(3,5),(5,7),(5,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 36 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(3,5),(5,7),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 37 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 38 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(3,5),(5,8),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 39 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(3,6),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 40 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(3,7),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=3, U=(6, 8), Y_U=(5, 8); fp = X_5_8/(X_1_3*X_1_4*X_1_8*X_3_7); cousins = 7 (all step-1 killable) |
| 41 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(3,8),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=3, U=(6, 8), Y_U=(5, 8); fp = X_5_8/(X_1_3*X_1_4*X_1_8*X_3_8); cousins = 7 (all step-1 killable) |
| 42 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(3,9),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 43 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,6),(4,9),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 44 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,6),(5,7),(5,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 45 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,6),(5,7),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 46 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,6),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 47 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,6),(5,8),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 48 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,6),(5,9),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{4,6}, k=3, U=(5, 9), Y_U=(4, 9); fp = X_4_9/(X_1_3*X_1_4*X_1_8*X_6_8); cousins = 7 (all step-1 killable) |
| 49 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,6),(6,8),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 50 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,6),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 51 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,7),(5,7),(5,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 52 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,7),(5,7),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 53 | CLUSTER | CROSSING | `{(1,3),(1,4),(1,8),(4,7),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 54 | CLUSTER | CROSSING | `{(1,3),(1,4),(2,8),(2,9),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 55 | CLUSTER | CROSSING | `{(1,3),(1,4),(2,9),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 56 | CLUSTER | CROSSING | `{(1,3),(1,4),(2,9),(4,6),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 57 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(2,9),(4,6),(6,8),(7,9)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=3, U=(6, 8), Y_U=(5, 8); fp = X_5_8/(X_1_3*X_1_4*X_2_9*X_7_9); cousins = 7 (all step-1 killable) |
| 58 | CLUSTER | CROSSING | `{(1,3),(1,4),(2,9),(4,6),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 59 | CLUSTER | CROSSING | `{(1,3),(1,4),(2,9),(4,7),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 60 | CLUSTER | CROSSING | `{(1,3),(1,4),(3,5),(3,9),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 61 | CLUSTER | CROSSING | `{(1,3),(1,4),(3,5),(4,9),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 62 | CLUSTER | CROSSING | `{(1,3),(1,4),(3,8),(3,9),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 63 | CLUSTER | CROSSING | `{(1,3),(1,4),(3,9),(4,6),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 64 | CLUSTER | CROSSING | `{(1,3),(1,4),(3,9),(4,6),(6,8),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 65 | CLUSTER | CROSSING | `{(1,3),(1,4),(3,9),(4,6),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 66 | CLUSTER | CROSSING | `{(1,3),(1,4),(3,9),(4,7),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 67 | CLUSTER | CROSSING | `{(1,3),(1,4),(4,6),(4,9),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 68 | NON_CLUSTER | CROSSING | `{(1,3),(1,4),(4,6),(4,9),(6,8),(7,9)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=3, U=(6, 8), Y_U=(5, 8); fp = X_5_8/(X_1_3*X_1_4*X_4_9*X_7_9); cousins = 7 (all step-1 killable) |
| 69 | CLUSTER | CROSSING | `{(1,3),(1,5),(1,8),(2,4),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 70 | NON_CLUSTER | CROSSING | `{(1,3),(1,5),(1,8),(3,5),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{7,9}, k=3, U=(1, 8), Y_U=(1, 7); fp = X_1_7/(X_1_3*X_1_5*X_3_5*X_4_6); cousins = 7 (all step-1 killable) |
| 71 | CLUSTER | CROSSING | `{(1,3),(1,5),(1,8),(3,5),(5,7),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 72 | CLUSTER | CROSSING | `{(1,3),(1,5),(1,8),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 73 | CLUSTER | CROSSING | `{(1,3),(1,5),(1,8),(3,6),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 74 | CLUSTER | CROSSING | `{(1,3),(1,5),(2,9),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 75 | CLUSTER | CROSSING | `{(1,3),(1,5),(3,5),(3,9),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 76 | CLUSTER | CROSSING | `{(1,3),(1,6),(1,7),(1,8),(2,4),(4,6)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 77 | CLUSTER | CROSSING | `{(1,3),(1,6),(1,7),(1,8),(3,5),(4,6)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 78 | CLUSTER | CROSSING | `{(1,3),(1,6),(1,7),(1,8),(3,5),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 79 | CLUSTER | CROSSING | `{(1,3),(1,6),(1,8),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 80 | CLUSTER | CROSSING | `{(1,3),(1,6),(2,4),(4,6),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 81 | CLUSTER | CROSSING | `{(1,3),(1,6),(2,9),(3,5),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 82 | CLUSTER | CROSSING | `{(1,3),(1,6),(3,5),(3,9),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 83 | CLUSTER | CROSSING | `{(1,3),(1,6),(3,5),(4,6),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 84 | CLUSTER | CROSSING | `{(1,3),(1,6),(3,5),(5,7),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 85 | CLUSTER | CROSSING | `{(1,3),(1,6),(3,5),(5,9),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 86 | CLUSTER | DOUBLE_POLE | `{(1,3),(1,7),(1,7),(1,8),(3,5),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 87 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(2,4),(3,5),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 88 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(2,4),(4,6),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 89 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(2,4),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 90 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(2,4),(4,7),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 91 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(2,5),(3,5),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 92 | NON_CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(2,6),(3,5),(5,7)}` | SINGLE-ORBIT CASCADE: Z_{1,3}, k=3, U=(2, 6), Y_U=(1, 6); fp = X_1_6/(X_1_7*X_1_8*X_3_5*X_5_7); cousins = 7 (all step-1 killable) |
| 93 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(2,7),(3,5),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 94 | NON_CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(2,9),(3,5),(5,7)}` | SINGLE-ORBIT CASCADE: Z_{4,6}, k=3, U=(5, 7), Y_U=(4, 7); fp = X_4_7/(X_1_3*X_1_7*X_1_8*X_2_9); cousins = 7 (all step-1 killable) |
| 95 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(3,5),(3,9),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 96 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(3,5),(4,6),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 97 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(3,5),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 98 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(3,5),(4,7),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 99 | NON_CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(3,5),(4,9),(5,7)}` | SINGLE-ORBIT CASCADE: Z_{2,4}, k=3, U=(3, 5), Y_U=(2, 5); fp = X_2_5/(X_1_7*X_1_8*X_4_9*X_5_7); cousins = 7 (all step-1 killable) |
| 100 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(3,5),(5,7),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 101 | NON_CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(3,5),(5,7),(6,9)}` | SINGLE-ORBIT CASCADE: Z_{2,4}, k=3, U=(3, 5), Y_U=(2, 5); fp = X_2_5/(X_1_7*X_1_8*X_5_7*X_6_9); cousins = 7 (all step-1 killable) |
| 102 | CLUSTER | CROSSING | `{(1,3),(1,7),(1,8),(3,6),(4,6),(5,7)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 103 | CLUSTER | CROSSING | `{(1,3),(1,8),(2,4),(2,9),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 104 | CLUSTER | CROSSING | `{(1,3),(1,8),(2,4),(3,5),(5,7),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 105 | NON_CLUSTER | CROSSING | `{(1,3),(1,8),(2,4),(3,5),(5,7),(7,9)}` | SINGLE-ORBIT CASCADE: Z_{6,8}, k=3, U=(7, 9), Y_U=(6, 9); fp = X_6_9/(X_1_3*X_1_8*X_2_4*X_3_5); cousins = 7 (all step-1 killable) |
| 106 | NON_CLUSTER | CROSSING | `{(1,3),(1,8),(2,4),(3,7),(4,6),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{5,7}, k=3, U=(6, 8), Y_U=(5, 8); fp = X_5_8/(X_1_3*X_1_8*X_2_4*X_3_7); cousins = 7 (all step-1 killable) |
| 107 | CLUSTER | CROSSING | `{(1,3),(1,8),(2,4),(3,9),(4,6),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 108 | CLUSTER | CROSSING | `{(1,3),(1,8),(2,4),(4,6),(5,7),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 109 | NON_CLUSTER | CROSSING | `{(1,3),(1,8),(2,4),(4,6),(5,9),(6,8)}` | SINGLE-ORBIT CASCADE: Z_{4,6}, k=3, U=(5, 9), Y_U=(4, 9); fp = X_4_9/(X_1_3*X_1_8*X_2_4*X_6_8); cousins = 7 (all step-1 killable) |
| 110 | CLUSTER | CROSSING | `{(1,3),(1,8),(2,4),(4,6),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 111 | CLUSTER | CROSSING | `{(1,3),(1,8),(2,4),(4,7),(5,7),(6,8)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 112 | CLUSTER | CROSSING | `{(1,3),(1,8),(3,5),(5,9),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |
| 113 | CLUSTER | CROSSING | `{(1,3),(1,8),(3,6),(4,6),(6,9),(7,9)}` | BLOCK-RULE: cluster matrix rank = 90 = full, so a_O = 0 once external columns are 0. |

## 3. Verification block

- **All reps non-local?**  triangulations = 0; crossings = 98; double-poles = 15; unclassified = 0.  YES ✓
- **Every NON_CLUSTER orbit has a verified single-orbit recipe?**  23/23 found; 0 timeouts, 0 no-recipe, 0 missing.  YES ✓
- **Cluster rank = full, nullity = 0?**  rank = 90, size = 90, nullity = 0.  YES ✓
- **Cluster + non-cluster = total?**  90 + 23 = 113;  total = 113.  YES ✓

### Overall verdict

> **n=9 LOCALITY: PROVEN.**

## 4. Per-orbit cascade table summary

Cascade outcomes restricted to NON_CLUSTER orbits (these are the orbits whose locality requires the single-orbit cascade to succeed):

| Outcome | Count |
|---|---:|
| recipe found | 23 |
| timed out    | 0 |
| exhausted, no recipe | 0 |
| missing record | 0 |

Cascade outcomes restricted to CLUSTER orbits (incidental — cluster orbits are killed jointly by the cluster matrix; their single-orbit cascade outcomes are reported here for completeness):

| Outcome | Count |
|---|---:|
| recipe found | 72 |
| timed out    | 14 |
| exhausted, no recipe | 4 |
| missing record | 0 |

All anomalies (timeouts + no-recipe) lie inside the cluster — so they are killed jointly by the block-rule, not by individual cascades.

## 5. Triangulation note (UNITARITY, not locality)

The 4 "untouched" Step-2 components in `../../step2_equate/flip_graph_n9/` are TRIANGULATION-only. They consist of fan-class triangulations (chord-length signature `(2, 2, 3, 3, 4, 4)`). Triangulations are LOCAL terms (= the Feynman-diagram coefficients of $A_n^{\text{tree}}$); they should NOT be killed by any mechanism in this folder. They survive Step-1 by design (see the local-survival lemma in `paper/step1 Kill technique and statistics/step-one-doesnt-kill-triangulations/`).

Equating the 49 triangulation $c$-values to a single common scalar $c$ is the *unitarity* claim. It is handled by the **d-subset argument** in `paper/Details missing in Hidden zero -> unitarity Rodina proof/details for Rodina D subset argument/` and `notes/11.` & `notes/18.`, independently of the cascade machinery here.

## 6. Conclusion

> **n=9 locality is fully proven from the cyclic 1-zero conditions alone**, by Step-1 (904 752 of 905 763 non-tri multisets killed directly) + 23 single-orbit cascades + the 90-orbit cluster matrix (rank 90 = full).

Unitarity at n=9 is established separately by the d-subset argument cited in §5 above.

