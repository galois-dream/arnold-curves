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
