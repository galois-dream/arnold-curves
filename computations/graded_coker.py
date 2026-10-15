"""Arnold's p.35 example (t^4, t^6 + t^7): graded cokernel and the first differential d_1.

Written 2026-10-15 14:50 CEST.  [Ours] reconstruction; [Arnold] p.35 names the Euler field 2x∂x+3y∂y and the
Hamiltonian field 2y∂x+3x^2∂y of the cusp x^3 = y^2.  Weights: t -> 1, x -> 4, y -> 6 (so the cusp's fields
appear as 4x∂x+6y∂y = 2*(2x∂x+3y∂y) and 2y∂x+3x^2∂y).  Weighted degree of t^m in slot i is m - a_i.
Run: .venv/bin/python computations/graded_coker.py
"""
from fractions import Fraction

import sympy as sp

from codim_check import K, Span, local_algebra, mul, series

t, x, y = sp.symbols("t x y")


def graded_coker(curve, lead):
    """dim of (F_d + TAf)/(F_{d+1} + TAf) for each weighted degree d, F_d = span of t^m e_i with m - a_i >= d."""
    fs = [series(c) for c in curve]
    N = len(fs)
    A = local_algebra(fs)
    gens = []
    for i in range(N):
        for r in A.rows.values():
            v = [Fraction(0)] * (N * K)
            v[i * K:(i + 1) * K] = r
            gens.append(v)
    dfs = [[m * f[m] for m in range(1, K)] + [Fraction(0)] for f in fs]
    for k in range(1, K):
        tk = [Fraction(0)] * K
        tk[k] = Fraction(1)
        gens.append([x_ for d in dfs for x_ in mul(tk, d)])
    mons = [(m - lead[i], i, m) for i in range(N) for m in range(1, K - 8)]

    def dim_with(dd):
        S = Span()
        for g in gens:
            S.add(g)
        for e, i, m in mons:
            if e >= dd:
                v = [Fraction(0)] * (N * K)
                v[i * K + m] = Fraction(1)
                S.add(v)
        return len(S.rows)

    out = {}
    for d in sorted({e for e, _, _ in mons}):
        c = dim_with(d) - dim_with(d + 1)
        if c:
            out[d] = c
    return out


def act(right, left, f):
    """Infinitesimal action of (right field right(t)·∂t, left field (left_x, left_y)) on f: right·f' - left∘f."""
    fx, fy = f
    return [sp.expand(right * sp.diff(fx, t) - left[0].subs({x: fx, y: fy}, simultaneous=True)),
            sp.expand(right * sp.diff(fy, t) - left[1].subs({x: fx, y: fy}, simultaneous=True))]


def toy_example():
    """[Ours] Appendix A toy: V = <v0 (deg 0), v1 (deg 1)>, W = <w0, w1, w2> (degrees 0, 1, 2);
    L v0 = w1 (raises degree by 1), L v1 = w1 + w2.  Graded cokernel pieces C_p = dim(F_p+LV) - dim(F_{p+1}+LV)."""
    LV = sp.Matrix([[0, 0], [1, 1], [0, 1]])          # columns L v0, L v1 in basis w0, w1, w2
    F = {p: sp.eye(3)[:, p:] for p in range(4)}         # F_p = span of basis vectors of degree >= p
    dim = lambda p: (F[p].row_join(LV)).rank() if p < 3 else LV.rank()
    C = [dim(p) - dim(p + 1) for p in range(3)]
    print("\nToy filtered map (Appendix A): rank L =", LV.rank(), "; graded cokernel C_0, C_1, C_2 =", C,
          "; L(v0 - v1) =", list(LV[:, 0] - LV[:, 1]), "(a degree-2 pivot found from a degree-0 vector: d_2)")


def eliminate(columns, col_deg, row_deg, pivots_out=None):
    """[Ours] Algorithm 3 of the note (added 2026-10-15 22:01 CEST): ordered elimination.
    columns: list of sparse vectors {row: coeff}; col_deg[k], row_deg[row]: degrees.  Columns are taken in
    DECREASING degree, so a column is only ever corrected by columns of higher or equal degree; each reduced
    column pivots on its lowest-degree row.  Returns ([(pivot row, r = deg(row) - deg(column))], unpivoted rows)."""
    lowest = lambda v: min(v, key=lambda q: (row_deg[q], q))
    pivots, gaps = ({} if pivots_out is None else pivots_out), []
    for k in sorted(range(len(columns)), key=lambda k: -col_deg[k]):
        v = {q: Fraction(c) for q, c in columns[k].items() if c}
        while v and lowest(v) in pivots:
            q = lowest(v)
            w = pivots[q]
            c = v[q] / w[q]
            for q2, c2 in w.items():
                v[q2] = v.get(q2, 0) - c * c2
                if v[q2] == 0:
                    del v[q2]
        if v:
            q = lowest(v)
            pivots[q] = v
            gaps.append((q, row_deg[q] - col_deg[k]))
    return gaps, [q for q in row_deg if q not in pivots]


def curve_moves(curve, lead):
    """Columns for eliminate(): right moves t^(j+1)d/dt (degree j) and left moves x^alpha d/dx_i
    (degree sum alpha_l d_l - d_i), as images mod t^K; rows (i, m) of degree m - d_i."""
    fs = [series(c) for c in curve]
    N = len(fs)
    row_deg = {(i, m): m - lead[i] for i in range(N) for m in range(1, K)}
    cols, degs = [], []
    for j in range(K - 1):  # right moves
        cols.append({(i, m + j): m * fs[i][m] for i in range(N) for m in range(1, K - j) if fs[i][m]})
        degs.append(j)
    def monos(start, used):  # exponent vectors alpha with sum alpha_l lead_l < K
        if start == N:
            yield ()
            return
        for e in range((K - 1 - used) // lead[start] + 1):
            for rest in monos(start + 1, used + e * lead[start]):
                yield (e,) + rest
    for alpha in monos(0, 0):
        if sum(alpha) == 0:
            continue
        g = [Fraction(0)] * K
        g[0] = Fraction(1)
        for l, e in enumerate(alpha):
            for _ in range(e):
                g = mul(g, fs[l])
        for i in range(N):
            cols.append({(i, m): c for m, c in enumerate(g) if c})
            degs.append(sum(e * d for e, d in zip(alpha, lead)) - lead[i])
    return cols, degs, row_deg


if __name__ == "__main__":
    for name, curve, lead in [("(t^4, t^6+t^7)", [4, (6, 7)], [4, 6]),
                              ("(t^4, t^6+t^9)", [4, (6, 9)], [4, 6]),
                              ("(t^4, t^6+t^7, t^9)", [4, (6, 7), 9], [4, 6, 9])]:
        g = graded_coker(curve, lead)
        print(f"{name:22} graded cokernel {g}  total {sum(g.values())}  "
              f"highest degree {max(g)}")
    f0, f = (t**4, t**6), (t**4, t**6 + t**7)
    print("\nFields that do nothing to the leading curve f0 = (t^4, t^6), but act on f = f0 + (0, t^7):")
    cases = [
        ("equation field (y^2 - x^3)∂y, degree 6", 0, (sp.Integer(0), y**2 - x**3)),
        ("Euler: t∂t paired with 4x∂x + 6y∂y, degree 0", t, (4 * x, 6 * y)),
        ("Hamiltonian: (1/2)t^3∂t paired with 2y∂x + 3x^2∂y, degree 2", t**3 / 2, (2 * y, 3 * x**2)),
    ]
    for label, right, left in cases:
        print(f"  {label}:\n     on f0: {act(right, left, f0)}   on f: {act(right, left, f)}")
    toy_example()
    print("\nAlgorithm 3 (ordered elimination) on the same curves: survivors by degree, and pivots by page")
    for name, curve, lead in [("(t^4, t^6+t^7)", [4, (6, 7)], [4, 6]), ("(t^4, t^6+t^9)", [4, (6, 9)], [4, 6]),
                              ("(t^4, t^6+t^7, t^9)", [4, (6, 7), 9], [4, 6, 9])]:
        gaps, free = eliminate(*curve_moves(curve, lead))
        surv = {}
        for (i, m) in free:
            if m < K - 8:
                surv[m - lead[i]] = surv.get(m - lead[i], 0) + 1
        rs = [r for (i, m), r in gaps if m < K - 8]  # pivots on rows below the truncation margin
        pages = {r: rs.count(r) for r in sorted(set(rs))}
        same = surv == graded_coker(curve, lead)
        print(f"  {name:22} survivors {dict(sorted(surv.items()))}  (= graded cokernel above: {same})  "
              f"pivots d_r by r: {pages}")
    toy_gaps, toy_free = eliminate([{"w1": 1}, {"w1": 1, "w2": 1}], [0, 1], {"w0": 0, "w1": 1, "w2": 2})
    print("  toy example: pivots (row, r) =", toy_gaps, "; survivors", toy_free)
