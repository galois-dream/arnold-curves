"""Small worked examples quoted in Sections 1-3 of the note.

Written 2026-10-15 14:51 CEST.  [Ours]
  (1) (t^2, t^3 + t^4) ~ (t^2, t^3) by one left change y -> y - x^2.
  (2) (t^3, t^4 + t^5) ~ (t^3, t^4) by explicit right/left changes (cf. [Arnold] p.30 Remark, E_{6,1,0} = (3,4)).
  (3) For f = (t^3, t^4): which t^{j+1} f' survive modulo A_f^2 (the b of Arnold's formula a N - b).
  (4) Sylvester duality for u=3, v=5 as quoted by [Arnold] p.33.
  (5) The rolling wheel of Figure 1, and (a, b) for (2,3), (3,4,5), (2,5), (3,4) [(2,5) added 2026-10-16 00:25 CEST]. Added 2026-10-15 21:21 CEST (corrected from a guessed 21:25).
Run: .venv/bin/python computations/note_examples.py
"""
import sympy as sp

from codim_check import a_b, local_algebra, series

t, T, S = sp.symbols("t T S")
ORD = 30
trunc = lambda e, v: sp.expand(sum(sp.expand(e).coeff(v, n) * v**n for n in range(ORD)))
orders = lambda e, v: sorted({m[0] for m in sp.Poly(e, v).monoms()}) if sp.expand(e) != 0 else []

# (1)
x, y = t**2, t**3 + t**4
print("(1) (t^2, t^3+t^4) after y -> y - x^2:", (x, sp.expand(y - x**2)))

# (2)
x, y = T**3, T**4 + T**5
a = sp.Rational(-1, 4)                                   # t = T(1 + aT): kills T^5 in y
x, y = trunc(x.subs(T, T * (1 + a * T)), T), trunc(y.subs(T, T * (1 + a * T)), T)
c4 = x.coeff(T, 4)
x = sp.expand(x - c4 * y)                                # left: kills T^4 in x
c5 = x.coeff(T, 5)
b = -c5 / 3                                              # T = S(1 + b S^2): kills S^5 in x
x, y = trunc(x.subs(T, S * (1 + b * S**2)), S), trunc(y.subs(T, S * (1 + b * S**2)), S)
c5y = y.coeff(S, 5)
print("(2) (t^3, t^4+t^5): remaining orders  x:", orders(x - S**3, S)[:5], " y:", orders(y - S**4, S)[:5],
      f" (y still has S^5 coefficient {c5y})")
cond = a_b([3, 4])[2]
print(f"    conductor of <3,4> = {cond}; remaining orders all >= {cond} and in the semigroup -> left-removable:",
      all(e >= cond for e in orders(x - S**3, S) + orders(y - S**4, S)))

# (3)
A = local_algebra([series(3), series(4)])
in_A = lambda m: m in A.rows                             # monomials t^m in A_f (A_f is monomial here)
print("(3) f = (t^3, t^4): A_f contains t^m for m in", sorted(A.rows)[:6], "...; gaps", [m for m in range(1, 8) if not in_A(m)])
for j in range(0, 5):
    comps = [(3, j + 3), (4, j + 4)]                     # t^{j+1} f' = (3 t^{j+3}, 4 t^{j+4})
    left = [f"{c}t^{e}" for c, e in comps if not in_A(e)]
    print(f"    t^{j + 1} f' = (3t^{j + 3}, 4t^{j + 4})  modulo A_f^2: {left or '0'}")
print("    a, b =", a_b([3, 4])[:2])

# (4)
u, v = 3, 5
K = (u - 1) * (v - 1)
rep = lambda c: any((c - n * v) >= 0 and (c - n * v) % u == 0 for n in range(c // v + 1))
row = [c for c in range(K)]
print(f"(4) u=3, v=5, K={K}: representable in [0,K-1]:", [c for c in row if rep(c)],
      " duality holds:", all(rep(c) != rep(K - 1 - c) for c in row))

# (5) Rolling wheel (Section 1, Figure 1): path (t - r sin t, 1 - r cos t).
r = sp.symbols("r", positive=True)
vel = (1 - r * sp.cos(t), r * sp.sin(t))
print("(5) wheel: velocity (1 - r cos t, r sin t) = 0 has solutions (r, t) =",
      sp.solve(list(vel), [r, t], dict=True), "(t mod 2pi) -> cusps only at r = 1")
print("    at r = 1: height", sp.series(1 - sp.cos(t), t, 0, 6).removeO(), "; position", sp.series(t - sp.sin(t), t, 0, 6).removeO())
print("    (2,3), (3,4,5), (2,5), (3,4): (a, b) =", a_b([2, 3])[:2], a_b([3, 4, 5])[:2], a_b([2, 5])[:2], a_b([3, 4])[:2],
      "-> codim N, 2N, 2N-1, 3N-2")
