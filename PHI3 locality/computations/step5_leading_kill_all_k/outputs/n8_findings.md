# Findings at n = 8

- self-check (k=1 only reproduces step-1 counts): **PASS**
- Phase 0 (existing step-1 non-tri survivors tested vs all (k,i)):
    - total tested: 100
    - killed by some higher-k zone: 100
    - **survives every (k,i)**: 0
    - elapsed: 0.2 s

- Phase 1 (full enumeration):
    - total multisets: 42504
    - triangulations: 132
    - non-tri killed by some (k,i): 42372
    - **non-tri surviving every (k,i)**: 0
    - elapsed: 5.1 s
