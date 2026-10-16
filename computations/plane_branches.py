r"""Plane branches: our counts against Hefez-Hernandes (2011) and the Ebey-Zariski elimination criteria.

Written 2026-10-16 15:24 CEST (corrected from a guessed 15:40).  [Lit] A. Hefez, M. E. Hernandes, "The analytic classification of plane
branches", Bull. London Math. Soc. 43 (2011) 289-298, p.291: Lambda = { ord_t(g(f) x' + h(f) y') + 1 } over
power series g, h in two variables; Gamma \ {0} lies in Lambda; Zariski invariant = min(Lambda \ Gamma) - m.
[Ours] (1) For coprime 3 <= m < n and F = mn - m - n: after x -> x + (m/n) y^(m-2) and the re-timing
s = x^(1/m), the curve (t^m, t^n + t^F) becomes (s^m, s^n + O(s^(F+1))).  Since s = t + higher, it suffices
that y(t) - s(t)^n = O(t^(F+1)); checked here symbolically.  (2) Arnold's re-timing number b (codim_check.a_b)
equals delta - 1 + #(Lambda \ Gamma) for every plane branch tried, delta = number of gaps.
(3) [Ours, added 2026-10-16 15:31 CEST] Zariski's theorem for (t^m, t^n + higher) without Puiseux theory:
keep x = t^m and remove the lowest term t^k of y by a left move in y (k in Gamma) or by a left move in x
followed by the re-timing s = x^(1/m) (k + m - n in Gamma).  If neither applies, the Euler vector
t f' - (m x, n y) = (0, (k - n) a t^k + ...) shows k + m is in Lambda \ Gamma.  Checked: the loop stops at k
exactly when Lambda \ Gamma is nonempty, then k + m = min(Lambda \ Gamma) (Zariski's invariant), and the
codimension drops below the monomial's by #(Lambda \ Gamma); otherwise the curve is the monomial.
Run: .venv/bin/python computations/plane_branches.py
"""
from fractions import Fraction
from itertools import combinations
from math import gcd

import sympy as sp

import codim_check
from codim_check import Span, a_b, codimension, local_algebra, mul, series

t = sp.symbols("t")


def frobenius_two_moves(m, n, upto=None):
    """Lowest order of y(t) - s(t)^n after the two moves (series cut at t^upto); the proof needs it > F."""
    F = m * n - m - n
    upto = upto or F + 2
    y = t**n + t**F
    x1 = t**m + sp.Rational(m, n) * y**(m - 2)
    sn = t**n * sp.series(sp.expand(x1 / t**m)**sp.Rational(n, m), t, 0, upto - n).removeO()
    rest = sp.expand(y - sn)
    return F, min(e for (e,) in sp.Poly(rest, t).monoms()) if rest != 0 else sp.oo


def lam(curve):
    """Lambda of a plane branch (x(t), y(t)), with the truncation certified by Algorithm 1."""
    a, b, K = codimension(curve)
    codim_check.K = K
    x, y = (series(c) for c in curve)
    A = local_algebra([x, y])
    dx, dy = ([e * f[e] if e < K else 0 for e in range(1, K)] + [0] for f in (x, y))
    one = [1] + [0] * (K - 1)
    D = Span()
    for g in [one] + list(A.rows.values()):
        D.add(mul(g, dx))
        D.add(mul(g, dy))
    gamma = set(A.rows) | {0}
    c = max(e for e in range(K) if e not in gamma) + 1
    L = {p + 1 for p in D.rows if p + 1 < c}  # every value >= c is in Gamma anyway
    codim_check.K = 48
    return a, b, gamma, L - gamma, c


def smul(u, v):
    w = [Fraction(0)] * len(u)
    for i, a in enumerate(u):
        if a:
            for j in range(len(u) - i):
                if v[j]:
                    w[i + j] += a * v[j]
    return w


def spow(u, r):
    """u**r for a truncated series with u[0] = 1 (binomial series)."""
    w, out = [Fraction(0)] + u[1:], [Fraction(1)] + [Fraction(0)] * (len(u) - 1)
    term, c = out[:], Fraction(1)
    for j in range(1, len(u)):
        c, term = c * (r - j + 1) / j, smul(term, w)
        if not any(term) or c == 0:
            break
        out = [o + c * x for o, x in zip(out, term)]
    return out


def compose(y, v):
    """y(s v(s)) with v[0] = 1."""
    out = [Fraction(0)] * len(y)
    for e, a in enumerate(y):
        if a:
            out[e:] = [o + a * x for o, x in zip(out[e:], spow(v, e))]
    return out


def rep(m, n, g):
    return next(((a, (g - a * m) // n) for a in range(g // m + 1) if (g - a * m) % n == 0), None)


def zariski_loop(m, n, extra):
    """Keep x = t^m; remove y's lowest term below the conductor while EC1 or EC2 applies.
    Returns (None, steps) if every term goes, else (k, steps) with t^k the first term that stays."""
    c = (m - 1) * (n - 1)
    y = [Fraction(0)] * (c + 1)
    y[n] = Fraction(1)
    for e, a in extra.items():
        y[e] += Fraction(a)
    tpow = lambda e: [Fraction(int(i == e)) for i in range(c + 1)]
    for step in range(4 * c):
        k = next((e for e in range(n + 1, c) if y[e]), None)
        if k is None:
            return None, step
        a = y[k]
        if rep(m, n, k):  # EC1: y -> y - a x^al y^be
            al, be = rep(m, n, k)
            y = [p - a * q for p, q in zip(y, smul(tpow(m * al), spow_int(y, be)))]
        elif rep(m, n, k + m - n):  # EC2: x -> x + (m a / n) x^al y^be, then s = x^(1/m)
            al, be = rep(m, n, k + m - n)
            q = smul(tpow(m * al), spow_int(y, be))
            u = spow([Fraction(int(i == 0)) + Fraction(m, n) * a * (q[i + m] if i + m <= c else 0)
                      for i in range(c + 1)], Fraction(1, m))  # s = t u(t)
            v = [Fraction(int(i == 0)) for i in range(c + 1)]
            for _ in range(c + 1):  # invert: t = s v(s)
                v = spow(compose(u, v), -1)
            y = compose(y, v)
        else:
            return k, step
    raise RuntimeError("loop did not finish")


def spow_int(y, b):
    out = [Fraction(int(i == 0)) for i in range(len(y))]
    for _ in range(b):
        out = smul(out, y)
    return out


def tjurina(curve):
    """tau = dim O_C / (f_X(f), f_Y(f)) O_C on truncated t-series.  The orders of an ideal J of O_C below K are
    exact; if e = min ord J, then e + Gamma lies in ord J, so every order >= e + c is in J; K >= e + c certifies."""
    m, ycomp = curve
    yt = sum(c * t**e for e, c in ({ycomp: 1} if isinstance(ycomp, int) else
             (ycomp if isinstance(ycomp, dict) else {e: 1 for e in ycomp})).items())
    X, Y = sp.symbols("X Y")
    f = sp.expand(sp.resultant(t**m - X, yt - Y, t))
    a, b, _ = codimension(curve)
    c = 2 * a
    K = codim_check.K = 2 * c + 2 * m + 4
    A = local_algebra([series(m), series(ycomp)])
    gens = []
    for g in (sp.diff(f, X), sp.diff(f, Y)):
        v = [Fraction(0)] * K
        for (e,), co in sp.Poly(sp.expand(g.subs({X: t**m, Y: yt})), t).terms():
            if e < K:
                v[e] += Fraction(int(co))
        gens.append(v)
    J = Span()
    for g in [[Fraction(1)] + [Fraction(0)] * (K - 1)] + list(A.rows.values()):
        for h in gens:
            J.add(mul(g, h))
    assert min(J.rows) + c <= K, "raise K"
    tau = len([e for e in set(A.rows) | {0} if e not in J.rows])
    codim_check.K = 48
    return a, b, tau


PLANE = ([[m, n] for m in range(2, 7) for n in range(m + 1, 14) if gcd(m, n) == 1]
         + [[3, (3 * k + 1, 3 * p + 2)] for k in range(1, 5) for p in range(k, 2 * k)]
         + [[3, (3 * p + 2, 3 * k + 1)] for p in range(1, 4) for k in range(p + 1, 2 * p + 1)]
         + [[4, (5, 7)], [4, (7, 9)], [4, (7, 10)], [4, (5, 6)], [5, (6, 7)], [5, (7, 8)], [5, (7, 9, 11)],
            [3, {7: 1, 8: 1, 11: 5}], [4, (9, 10, 11)]]
         + [[4, (6, 2 * k + 1)] for k in range(3, 7)] + [[4, (6, 7, 9)], [6, (9, 10)], [4, (10, 11)]])

if __name__ == "__main__":
    print("(1) Two moves remove the Frobenius term of (t^m, t^n + t^F): lowest leftover order vs F")
    bad = 0
    for m in range(3, 7):
        for n in range(m + 1, 16):
            if gcd(m, n) == 1:
                F, low = frobenius_two_moves(m, n)
                bad += not low > F
    print(f"    all coprime 3 <= m <= 6 < n <= 15: {'every leftover is past F' if not bad else f'{bad} FAILURES'}")
    for m, n in ((3, 4), (3, 7), (4, 5), (5, 7)):
        F = m * n - m - n
        _, low = frobenius_two_moves(m, n, upto=2 * F + 2)
        print(f"    (t^{m}, t^{n} + t^{F}): leftover starts at order {low} (2F - n = {2 * F - n}), conductor {F + 1}")

    print("(2) b = delta - 1 + #(Lambda \\ Gamma) on plane branches")
    bad = 0
    for curve in PLANE:
        a, b, gamma, extra, c = lam(curve)
        ok = b == a - 1 + len(extra)
        bad += not ok
        if not ok or curve in ([3, 7], [3, (7, 8)], [4, (6, 7)], [3, {7: 1, 8: 1, 11: 5}]):
            zar = f", Zariski invariant {min(extra) - curve[0]}" if extra else ""
            print(f"    {str(curve):24} delta={a:<3} b={b:<3} Lambda\\Gamma={sorted(extra)}{zar}"
                  + ("" if ok else "   <-- MISMATCH"))
    print(f"    {len(PLANE)} plane branches, {bad} mismatches")

    print("(3) Zariski's theorem by elimination: where the loop stops vs Lambda \\ Gamma and the codimension")
    bad = total = 0
    for m, n in ((3, 4), (3, 5), (3, 7), (3, 8), (4, 5), (4, 7), (5, 6), (5, 7)):
        c = (m - 1) * (n - 1)
        extras = [{k: 1} for k in range(n + 1, c)] + [{k1: 1, k2: 2} for k1, k2 in combinations(range(n + 1, c), 2)
                                                        if m == 3 or k1 % 3 == 0]
        for ex in extras:
            k, steps = zariski_loop(m, n, ex)
            curve = [m, {n: 1, **ex}]
            a, b, gamma, extra, _ = lam(curve)
            ok = (extra == set()) if k is None else (min(extra) == k + m and b - (a - 1) == len(extra))
            total += 1
            bad += not ok
            if not ok or (m, n, tuple(ex)) in ((3, 7, (8,)), (3, 7, (11,)), (4, 7, (9,)), (4, 7, (13,))):
                print(f"    (t^{m}, t^{n} + {ex}): " + ("monomial" if k is None else f"stops at t^{k}, Zariski invariant {k}")
                      + f" after {steps} moves; Lambda\\Gamma={sorted(extra)}" + ("" if ok else "   <-- MISMATCH"))
    print(f"    {total} curves, {bad} mismatches")

    print("(4) Codimension at N = 3 against the Tjurina number tau (Zariski 1966, Thm 1 and Note)")
    bad = 0
    TAU = [[3, 4], [3, 7], [3, (7, 8)], [3, (7, 10)], [3, (10, 11)], [3, (10, 14)], [4, (5, 7)], [4, (7, 9)],
           [4, (7, 10)], [4, (6, 7)], [5, (6, 7)], [5, (7, 8)], [3, {7: 1, 8: 1, 11: 5}], [4, (6, 7, 9)]]
    for curve in TAU:
        a, b, tau = tjurina(curve)
        ok = 3 * a - b == tau + 1
        bad += not ok
        if not ok or curve in ([3, 7], [3, (7, 8)], [4, (6, 7)]):
            print(f"    {str(curve):12} codim at N=3: {3 * a - b}, tau = {tau}" + ("" if ok else "   <-- MISMATCH"))
    print(f"    {len(TAU)} plane branches, {bad} mismatches")
