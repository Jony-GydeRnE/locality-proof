# n=9 cluster: chord-length-spectrum reduction

Coarser than D_9.  Groups multisets by their sorted multiset of chord lengths.

- Z_9 orbits in cluster: **90**
- D_9 groups: **55**
- **Length-spectrum shapes: 8**
- Total multisets covered: 804

## The rule

A chord {a, b} of the 9-gon has LENGTH = min(|a-b|, 9-|a-b|), so length ∈ {2, 3, 4}. The LENGTH SPECTRUM of a 6-chord multiset is the sorted multiset of chord lengths. Two multisets are shape-equivalent iff they have the same length spectrum.

Length is preserved by every rotation and reflection of the 9-gon, so this is coarser than D_9: every D_9 orbit lives inside one length-spectrum shape, but a single shape can contain multiple D_9 orbits at different polygon positions.

## All shapes

| Shape # | Length spectrum | # D_9 groups | # Z_9 orbits | # multisets |
|---:|---|---:|---:|---:|
| 1 | 2^6 | 5 | 7 | 57 |
| 2 | 2^5 · 3 | 14 | 26 | 234 |
| 3 | 2^5 · 4 | 4 | 6 | 54 |
| 4 | 2^4 · 3^2 | 17 | 28 | 252 |
| 5 | 2^4 · 3 · 4 | 11 | 16 | 144 |
| 6 | 2^3 · 3^3 | 2 | 4 | 36 |
| 7 | 2^3 · 3^2 · 4 | 1 | 2 | 18 |
| 8 | 2^3 · 3 · 4^2 | 1 | 1 | 9 |

## Shape members (which D_9 groups land in each shape)

### Shape 1 — spectrum 2^6

| D_9 grp rep | Kind | Z_9 orbits | Multiset |
|---:|---|---|---|
| 4 | reflection_pair | 4, 12 | `{(1,3), (1,3), (1,8), (2,4), (4,6), (6,8)}` |
| 6 | reflection_pair | 6, 9 | `{(1,3), (1,3), (1,8), (3,5), (4,6), (6,8)}` |
| 103 | palindromic | 103 | `{(1,3), (1,8), (2,4), (2,9), (4,6), (6,8)}` |
| 104 | palindromic_leak | 104 | `{(1,3), (1,8), (2,4), (3,5), (5,7), (6,8)}` |
| 108 ★ | palindromic | 108 | `{(1,3), (1,8), (2,4), (4,6), (5,7), (7,9)}` |

### Shape 2 — spectrum 2^5 · 3

| D_9 grp rep | Kind | Z_9 orbits | Multiset |
|---:|---|---|---|
| 1 | reflection_pair | 1, 13 | `{(1,3), (1,3), (1,4), (1,8), (4,6), (6,8)}` |
| 2 | reflection_pair | 2, 5 | `{(1,3), (1,3), (1,7), (1,8), (3,5), (5,7)}` |
| 3 | reflection_pair | 3, 11 | `{(1,3), (1,3), (1,7), (3,5), (5,7), (7,9)}` |
| 7 | reflection_pair | 7, 10 | `{(1,3), (1,3), (1,8), (3,5), (5,7), (5,8)}` |
| 28 | palindromic | 28 | `{(1,3), (1,4), (1,8), (2,4), (4,6), (6,8)}` |
| 34 | reflection_pair | 34, 100 | `{(1,3), (1,4), (1,8), (3,5), (4,6), (6,8)}` |
| 36 | reflection_pair | 36, 97 | `{(1,3), (1,4), (1,8), (3,5), (5,7), (6,8)}` |
| 37 | reflection_pair | 37, 89 | `{(1,3), (1,4), (1,8), (3,5), (5,7), (7,9)}` |
| 45 | reflection_pair | 45, 96 | `{(1,3), (1,4), (1,8), (4,6), (5,7), (6,8)}` |
| 46 ★ | reflection_pair | 46, 88 | `{(1,3), (1,4), (1,8), (4,6), (5,7), (7,9)}` |
| 49 | reflection_pair | 49, 87 | `{(1,3), (1,4), (1,8), (4,6), (6,8), (7,9)}` |
| 55 | reflection_pair | 55, 107 | `{(1,3), (1,4), (2,9), (3,5), (5,7), (7,9)}` |
| 56 | reflection_pair | 56, 111 | `{(1,3), (1,4), (2,9), (4,6), (5,7), (7,9)}` |
| 110 | palindromic_leak | 110 | `{(1,3), (1,8), (2,4), (4,6), (6,9), (7,9)}` |

### Shape 3 — spectrum 2^5 · 4

| D_9 grp rep | Kind | Z_9 orbits | Multiset |
|---:|---|---|---|
| 69 | reflection_pair | 69, 79 | `{(1,3), (1,5), (1,8), (2,4), (4,6), (6,8)}` |
| 71 | reflection_pair | 71, 72 | `{(1,3), (1,5), (1,8), (3,5), (5,7), (6,8)}` |
| 74 | palindromic_leak | 74 | `{(1,3), (1,5), (2,9), (3,5), (5,7), (7,9)}` |
| 81 | palindromic_leak | 81 | `{(1,3), (1,6), (2,9), (3,5), (5,7), (7,9)}` |

### Shape 4 — spectrum 2^4 · 3^2

| D_9 grp rep | Kind | Z_9 orbits | Multiset |
|---:|---|---|---|
| 14 | reflection_pair | 14, 86 | `{(1,3), (1,4), (1,4), (1,8), (4,6), (6,8)}` |
| 21 | reflection_pair | 21, 23 | `{(1,3), (1,4), (1,7), (1,8), (3,5), (5,7)}` |
| 22 ★ | palindromic | 22 | `{(1,3), (1,4), (1,7), (1,8), (4,6), (5,7)}` |
| 24 | reflection_pair | 24, 98 | `{(1,3), (1,4), (1,7), (3,5), (5,7), (7,9)}` |
| 25 | reflection_pair | 25, 90 | `{(1,3), (1,4), (1,7), (4,6), (5,7), (7,9)}` |
| 35 | reflection_pair | 35, 91 | `{(1,3), (1,4), (1,8), (3,5), (5,7), (5,8)}` |
| 38 | reflection_pair | 38, 47 | `{(1,3), (1,4), (1,8), (3,5), (5,8), (6,8)}` |
| 39 | palindromic | 39 | `{(1,3), (1,4), (1,8), (3,6), (4,6), (6,8)}` |
| 42 | reflection_pair | 42, 60 | `{(1,3), (1,4), (1,8), (3,9), (4,6), (6,8)}` |
| 44 | reflection_pair | 44, 102 | `{(1,3), (1,4), (1,8), (4,6), (5,7), (5,8)}` |
| 50 | reflection_pair | 50, 58 | `{(1,3), (1,4), (1,8), (4,6), (6,9), (7,9)}` |
| 52 | reflection_pair | 52, 59 | `{(1,3), (1,4), (1,8), (4,7), (5,7), (6,8)}` |
| 53 | palindromic | 53 | `{(1,3), (1,4), (1,8), (4,7), (5,7), (7,9)}` |
| 54 | palindromic | 54 | `{(1,3), (1,4), (2,8), (2,9), (4,6), (6,8)}` |
| 63 | reflection_pair | 63, 64 | `{(1,3), (1,4), (3,9), (4,6), (5,7), (7,9)}` |
| 95 | palindromic_leak | 95 | `{(1,3), (1,7), (1,8), (3,5), (3,9), (5,7)}` |
| 113 | palindromic_leak | 113 | `{(1,3), (1,8), (3,6), (4,6), (6,9), (7,9)}` |

### Shape 5 — spectrum 2^4 · 3 · 4

| D_9 grp rep | Kind | Z_9 orbits | Multiset |
|---:|---|---|---|
| 16 | reflection_pair | 16, 77 | `{(1,3), (1,4), (1,5), (1,8), (5,7), (6,8)}` |
| 17 | reflection_pair | 17, 76 | `{(1,3), (1,4), (1,5), (1,8), (5,7), (7,9)}` |
| 43 | reflection_pair | 43, 93 | `{(1,3), (1,4), (1,8), (4,6), (4,9), (6,8)}` |
| 61 | reflection_pair | 61, 84 | `{(1,3), (1,4), (3,5), (4,9), (5,7), (7,9)}` |
| 67 | reflection_pair | 67, 83 | `{(1,3), (1,4), (4,6), (4,9), (5,7), (7,9)}` |
| 73 | palindromic_leak | 73 | `{(1,3), (1,5), (1,8), (3,6), (4,6), (6,8)}` |
| 75 | palindromic_leak | 75 | `{(1,3), (1,5), (3,5), (3,9), (5,7), (7,9)}` |
| 78 | palindromic_leak | 78 | `{(1,3), (1,6), (1,7), (1,8), (3,5), (5,7)}` |
| 80 | palindromic_leak | 80 | `{(1,3), (1,6), (2,4), (4,6), (6,9), (7,9)}` |
| 82 | palindromic_leak | 82 | `{(1,3), (1,6), (3,5), (3,9), (5,7), (7,9)}` |
| 112 | palindromic_leak | 112 | `{(1,3), (1,8), (3,5), (5,9), (6,9), (7,9)}` |

### Shape 6 — spectrum 2^3 · 3^3

| D_9 grp rep | Kind | Z_9 orbits | Multiset |
|---:|---|---|---|
| 27 | reflection_pair | 27, 65 | `{(1,3), (1,4), (1,7), (4,6), (6,9), (7,9)}` |
| 51 | reflection_pair | 51, 66 | `{(1,3), (1,4), (1,8), (4,7), (5,7), (5,8)}` |

### Shape 7 — spectrum 2^3 · 3^2 · 4

| D_9 grp rep | Kind | Z_9 orbits | Multiset |
|---:|---|---|---|
| 19 | reflection_pair | 19, 62 | `{(1,3), (1,4), (1,5), (3,9), (5,7), (7,9)}` |

### Shape 8 — spectrum 2^3 · 3 · 4^2

| D_9 grp rep | Kind | Z_9 orbits | Multiset |
|---:|---|---|---|
| 85 | palindromic_leak | 85 | `{(1,3), (1,6), (3,5), (5,9), (6,9), (7,9)}` |

★ = the D_9 group contains a cascade-failure seed orbit (22, 46, 88, or 108).
