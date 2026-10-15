"""Gibson–Hobbs (1993) Tables 1–2 checked against Arnold's codimension aN - b at N = 3.

Written 2026-10-15 18:35 CEST (corrected from a guessed 18:52).  [Lit] Gibson–Hobbs, Math. Proc. Camb. Phil. Soc. 113 (1993) 297–310:
Table 1 (planar, p.298) and Table 2 (spatial, p.299) give normal forms of the A-simple germs C -> C^3,
with value-semigroup generators S(f) and A_e-codimension c(f).  By their p.303 ("the A_e-codimension of f
differs from the codimension of the A^(k)-orbit by 1"), the orbit codimension is c(f) + 1.
[Ours] For each normal form, over a range of indices: does 3a - b (computations/codim_check.py) equal
c(f) + 1, and do the semigroup generators match S(f)?  Then: which rows share (semigroup, codimension)?
Run: .venv/bin/python computations/gibson_hobbs_check.py
"""
from collections import defaultdict

from codim_check import a_b, local_algebra, series

# (label, curve, S(f), c(f)) per index; binomial (p, q) means t^p + t^q; planar rows omit the zero component.
ROWS = []
for k in range(1, 6):  # Table 1
    ROWS.append(("T1 (t^2,t^(2k+1),0)", [2, 2 * k + 1], {2, 2 * k + 1}, 2 * k))
for k in range(1, 5):
    ROWS.append(("T1 (t^3,t^(3k+1),0)", [3, 3 * k + 1], {3, 3 * k + 1}, 6 * k))
    ROWS.append(("T1 (t^3,t^(3k+2),0)", [3, 3 * k + 2], {3, 3 * k + 2}, 6 * k + 2))
    for p in range(k, 2 * k):  # k <= p < 2k
        ROWS.append(("T1 (t^3,t^(3k+1)+t^(3p+2),0)", [3, (3 * k + 1, 3 * p + 2)], {3, 3 * k + 1}, 4 * k + p + 1))
for p in range(1, 4):
    for k in range(p + 1, 2 * p + 1):  # p < k <= 2p
        ROWS.append(("T1 (t^3,t^(3k+1)+t^(3p+2),0)", [3, (3 * p + 2, 3 * k + 1)], {3, 3 * p + 2}, 4 * p + k + 2))
ROWS += [("T1 (t^4,t^5,0)", [4, 5], {4, 5}, 12), ("T1 (t^4,t^5+t^7,0)", [4, (5, 7)], {4, 5}, 11),
         ("T1 (t^4,t^7,0)", [4, 7], {4, 7}, 18), ("T1 (t^4,t^7+t^9,0)", [4, (7, 9)], {4, 7}, 16),
         ("T1 (t^4,t^7+t^10,0)", [4, (7, 10)], {4, 7}, 17)]
for k in range(3, 6):
    ROWS.append(("T1 (t^4,t^6+t^(2k+1),0)", [4, (6, 2 * k + 1)], {4, 6, 2 * k + 7}, 2 * k + 8))
for k in range(1, 5):  # Table 2
    for n in range(k, 2 * k):  # k <= n < 2k
        ROWS.append(("T2 (t^3,t^(3k+1),t^(3n+2))", [3, 3 * k + 1, 3 * n + 2], {3, 3 * k + 1, 3 * n + 2}, 2 * k + 2 * n + 1))
        for p in range(k, n):  # k <= p < n < 2k
            ROWS.append(("T2 (t^3,t^(3k+1)+t^(3p+2),t^(3n+2))", [3, (3 * k + 1, 3 * p + 2), 3 * n + 2],
                         {3, 3 * k + 1, 3 * n + 2}, 2 * k + p + n + 1))
for n in range(1, 4):
    for k in range(n + 1, 2 * n + 1):  # n < k <= 2n
        ROWS.append(("T2 (t^3,t^(3n+2),t^(3k+1))", [3, 3 * n + 2, 3 * k + 1], {3, 3 * n + 2, 3 * k + 1}, 2 * k + 2 * n + 1))
        for p in range(n + 1, k):  # n < p < k
            ROWS.append(("T2 (t^3,t^(3n+2)+t^(3p+1),t^(3k+1))", [3, (3 * n + 2, 3 * p + 1), 3 * k + 1],
                         {3, 3 * n + 2, 3 * k + 1}, 2 * n + p + k + 1))
ROWS += [("T2 (t^4,t^5,t^6)", [4, 5, 6], {4, 5, 6}, 8), ("T2 (t^4,t^5,t^7)", [4, 5, 7], {4, 5, 7}, 9),
         ("T2 (t^4,t^5,t^11)", [4, 5, 11], {4, 5, 11}, 11), ("T2 (t^4,t^5+t^7,t^11)", [4, (5, 7), 11], {4, 5, 11}, 10),
         ("T2 (t^4,t^7,t^9+t^10)", [4, 7, (9, 10)], {4, 7, 9}, 12), ("T2 (t^4,t^7,t^9)", [4, 7, 9], {4, 7, 9}, 13),
         ("T2 (t^4,t^7,t^10)", [4, 7, 10], {4, 7, 10}, 14), ("T2 (t^4,t^7,t^13)", [4, 7, 13], {4, 7, 13}, 15),
         ("T2 (t^4,t^7,t^17)", [4, 7, 17], {4, 7, 17}, 17), ("T2 (t^4,t^7+t^9,t^10)", [4, (7, 9), 10], {4, 7, 10}, 13),
         ("T2 (t^4,t^7+t^9,t^13)", [4, (7, 9), 13], {4, 7, 13}, 14), ("T2 (t^4,t^7+t^9,t^17)", [4, (7, 9), 17], {4, 7, 17}, 15),
         ("T2 (t^4,t^7+t^10,t^17)", [4, (7, 10), 17], {4, 7, 17}, 16)]
for k in range(3, 7):
    ROWS.append(("T2 (t^4,t^6,t^(2k+1))", [4, 6, 2 * k + 1], {4, 6, 2 * k + 1}, 2 * k + 4))
    if k >= 4:
        ROWS.append(("T2 (t^4,t^6+t^(2k-1),t^(2k+1))", [4, (6, 2 * k - 1), 2 * k + 1], {4, 6, 2 * k + 1}, 2 * k + 3))
    if k >= 5:
        ROWS.append(("T2 (t^4,t^6+t^(2k-3),t^(2k+1))", [4, (6, 2 * k - 3), 2 * k + 1], {4, 6, 2 * k + 1, 2 * k + 3}, 2 * k + 2))
for k in range(7, 10):
    ROWS.append(("T2 (t^4,t^6+t^(2k-7),t^(2k+1))", [4, (6, 2 * k - 7), 2 * k + 1], {4, 6, 2 * k - 1, 2 * k + 1}, 2 * k - 1))


def generators(curve):
    A = local_algebra([series(c) for c in curve])
    S = sorted(A.rows)  # nonzero semigroup elements below the truncation
    return {s for s in S if not any(s - u in A.rows for u in S if u < s)}, frozenset(S)


def show(curve):
    return "(" + ",".join(str(c) if isinstance(c, int) else "(" + ",".join(map(str, c)) + ")" for c in curve) + ")"


if __name__ == "__main__":
    bad, seen = 0, defaultdict(list)
    for label, curve, S_gh, c_gh in ROWS:
        a, b, _ = a_b(curve)
        gens, S = generators(curve)
        ok = 3 * a - b == c_gh + 1 and gens == S_gh
        bad += not ok
        seen[(S, 3 * a - b)].append(show(curve))
        if not ok:
            print(f"MISMATCH {label} {show(curve)}: 3a-b={3 * a - b} vs c+1={c_gh + 1}; gens {sorted(gens)} vs {sorted(S_gh)}")
    print(f"{len(ROWS)} normal forms checked; {bad} disagree with 3a-b = c(f)+1 or with S(f).")
    print("Normal forms sharing semigroup and codimension (not separated by these invariants):")
    for (S, cod), cs in seen.items():
        if len(cs) > 1:
            print(f"  codim {cod}: {', '.join(cs)}")
