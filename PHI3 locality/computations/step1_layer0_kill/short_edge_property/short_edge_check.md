# Does every step-1 non-tri survivor have a short isolated chord?

**Property checked:** ∃ chord c ∈ multiset M with length 2 such that c does NOT cross any other chord in M.  (Length 2 is the shortest chord on an n-gon.  Sharing an endpoint with another chord does not count as crossing.)

## Per-n results

| n | non-tri survivors | with property | without | universal? |
|---:|---:|---:|---:|:---:|
| 7 | 7 | 7 | 0 | ✅ YES |
| 8 | 100 | 100 | 0 | ✅ YES |
| 9 | 1011 | 1008 | 3 | ❌ NO |

## n = 9 — counterexamples

3 multisets have NO length-2 chord that avoids crossings.  First 3 shown:

- `{(1,3), (1,8), (2,4), (4,6), (5,7), (7,9)}`
- `{(1,3), (2,9), (3,5), (4,6), (6,8), (7,9)}`
- `{(1,8), (2,4), (2,9), (3,5), (5,7), (6,8)}`

