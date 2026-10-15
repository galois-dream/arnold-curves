"""The class starting (t^3, t^10): Castro–Montgomery's 10 germs vs Arnold's E_18 list.

Written 2026-10-15 14:48 CEST.
[Lit] Castro–Montgomery, arXiv 0912.3020v1, p.12 (Thm 3.27 example): lists 10 germs, "precisely 10 orbits".
[Arnold] p.30, Thm 4 + Remark: E_{6k,p,i}, codim (2k+p)N - (2k+2p-2-i); boundary forms; E_{6k,k,k-1} = (3,3k+1).
[Ours] (1) codim aN-b of each germ vs Arnold's formula (k=3);
       (2) explicit change of variables showing (t^3, t^10 + t^17) ~ (t^3, t^10).
Run: .venv/bin/python computations/cm_3_10_orbits.py
"""
import sympy as sp

from codim_check import a_b

# germ -> (p, i) in Arnold's E_{18,p,i}; p = 3 = k means planar (third component absent)
GERMS = {
    "(3,10,11)": ([3, 10, 11], 0, 0),
    "(3,10+11,14)": ([3, (10, 11), 14], 1, 0),
    "(3,10,14)": ([3, 10, 14], 1, 1),
    "(3,10+11,17)": ([3, (10, 11), 17], 2, 0),
    "(3,10+14,17)": ([3, (10, 14), 17], 2, 1),
    "(3,10,17)": ([3, 10, 17], 2, 2),
    "(3,10+11)": ([3, (10, 11)], 3, 0),
    "(3,10+14)": ([3, (10, 14)], 3, 1),
    "(3,10+17)": ([3, (10, 17)], 3, 2),
    "(3,10)": ([3, 10], 3, 2),  # Arnold's Remark: E_{6k,k,k-1} = (3,3k+1)
}
k = 3
print("(1) codimension aN-b vs Arnold's E_18 formula (2k+p)N-(2k+2p-2-i):")
for name, (curve, p, i) in GERMS.items():
    a, b, c = a_b(curve)
    A, B = 2 * k + p, 2 * k + 2 * p - 2 - i
    print(f"  {name:14} computed {a}N-{b:<3} Arnold E_18,{p},{i}: {A}N-{B:<3} {'ok' if (a, b) == (A, B) else 'MISMATCH'}")

print("\n(2) explicit equivalence (t^3, t^10 + t^17) ~ (t^3, t^10):")
T, S = sp.symbols("T S")
ORD = 40
x, y = T**3, T**10 + T**17
a = sp.Rational(-1, 10)                    # right change t = T(1 + a T^7) kills T^17 in y
x = sp.expand(x.subs(T, T * (1 + a * T**7)))
y = sp.expand(y.subs(T, T * (1 + a * T**7)))
x = sp.expand(x - 3 * a * y)               # left change x -> x - 3a y kills T^10 in x
c17 = sp.Poly(x, T).coeff_monomial(T**17)
b = -c17 / 3                               # right change T = S(1 + b S^14) kills S^17 in x
x = sp.expand(x.subs(T, S * (1 + b * S**14)))
y = sp.expand(y.subs(T, S * (1 + b * S**14)))
trunc = lambda e: sum(e.coeff(S, n) * S**n for n in range(ORD))
rx, ry = sp.expand(trunc(x) - S**3), sp.expand(trunc(y) - S**10)
ex = sorted({m[0] for m in sp.Poly(rx, S).monoms()}) if rx != 0 else []
ey = sorted({m[0] for m in sp.Poly(ry, S).monoms()}) if ry != 0 else []
print(f"  after 2 right + 1 left changes: x = S^3 + (terms of orders {ex[:4]}...), y = S^10 + (orders {ey[:4]}...)")
cond = a_b([3, 10])[2]
ok = all(e >= cond for e in ex + ey)
print(f"  conductor of <3,10> = {cond}; all remaining orders >= conductor: {ok}")
print("  => remaining terms lie in A_f and are removed by left changes ([Arnold] p.34): EQUIVALENT" if ok else "  => NOT shown")
