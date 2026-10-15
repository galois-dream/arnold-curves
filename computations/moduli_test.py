"""Moduli inside a symbol class, by a dimension count; saturation numbers; an index-range check on Arnold Thm 4.

Written 2026-10-15 ~15:08 CEST (approx.; committed 15:16).  [Ours]
Symbol (Castro–Montgomery 2012, Def. 3.28): [m,n,p] = curves equivalent to (t^m,t^n,t^p,0,...) + O(t^{p+1});
[m,(n,p)] likewise with t^n + t^p.  Such a class S is the preimage of ONE orbit of the p-jet group in J^p,
so codim S = codim of that jet orbit = N*p - rank of the tangent space TA j projected to orders <= p.
If even a generic germ in S has an orbit of larger codimension, finitely many orbits cannot fill S
(S is irreducible: preimage of an orbit of a connected group), so S contains moduli.  The first jet order r
at which the generic r-jet orbit is thinner than S is where the modulus appears.
Run: .venv/bin/python computations/moduli_test.py [A|B|C]   (A takes ~2 min)
"""
import random
from fractions import Fraction

from codim_check import K, Span, a_b, local_algebra, mul, series


def tangent_gens(fs):
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
        gens.append([x for d in dfs for x in mul(tk, d)])
    return gens


def jet_codim(fs, r):
    """codim of the orbit of the r-jet of f in the space of r-jets (orders 1..r in each of N slots)."""
    N = len(fs)
    S = Span()
    for g in tangent_gens(fs):
        S.add([g[i * K + m] for i in range(N) for m in range(1, r + 1)])
    return N * r - len(S.rows)


def generic(symbol, N, p, top, seed):
    """symbol completed by random rational terms of orders p+1..top in every slot."""
    rnd = random.Random(seed)
    fs = []
    for i in range(N):
        terms = {}
        if i < len(symbol):
            c = symbol[i]
            for e in ((c,) if isinstance(c, int) else c):
                terms[e] = 1
        for e in range(p + 1, top + 1):
            terms[e] = terms.get(e, 0) + Fraction(rnd.randint(-9, 9), rnd.randint(1, 9))
        fs.append(series(terms))
    return fs


def moduli(symbol, N, p, top=30, seeds=(1, 2)):
    sym = [series(c) for c in symbol] + [series({})] * (N - len(symbol))
    cS = jet_codim(sym, p)
    out = []
    for sd in seeds:
        g = generic(symbol, N, p, top, sd)
        prof = {r: jet_codim(g, r) for r in range(p, 36)}
        first = next((r for r in sorted(prof) if prof[r] > cS), None)
        out.append((prof[35], first))
    return cS, out


if __name__ == "__main__":
  import sys
  parts = sys.argv[1] if len(sys.argv) > 1 else "ABC"
  if "A" in parts:
    print("A. Moduli test: codim of symbol class S vs codim of a generic orbit inside S (N = 3 and 4).")
    tests = [  # (label, symbol, p, CM2012 status)
        ("[3,7]", [3, 7], 7, "R^2CT: tower simple (CM Thm 3.26)"),
        ("[4,6,7]", [4, 6, 7], 7, "RCL: tower simple"),
        ("[4,6,9]", [4, 6, 9], 9, "RCRC: tower simple"),
        ("[4,(6,9)]", [4, (6, 9)], 9, "RCRC: tower simple"),
        ("[4,10,11]", [4, 10, 11], 11, "R^2CL: moduli claimed"),
        ("[4,10,13]", [4, 10, 13], 13, "R^2CRC: moduli claimed"),
        ("[4,(10,13)]", [4, (10, 13)], 13, "R^2CRC: moduli claimed"),
    ]
    for label, sym, p, status in tests:
        for N in (3, 4):
            cS, res = moduli(sym, N, p)
            gen = ", ".join(f"generic {c} (moduli >= {c - cS}; first at jet order {r})" for c, r in res)
            print(f"  {label:12} N={N}: codim S = {cS:3};  {gen}   [{status}]")

  if "B" in parts:
    print("\nB. Saturation numbers N(omega, inf) = number of curve orbits with symbol [3, 3k+1] / [3, 3k+2],")
    print("   from Arnold's Thm 4 list (p.30) read with i <= p, p <= k, and (3,(3k+1,6k-1)) ~ (3,3k+1):")
    for k in range(1, 5):
        for name, n, base in (("[3,3k+1]", 3 * k + 1, 3 * k + 2), ("[3,3k+2]", 3 * k + 2, 3 * k + 4)):
            orbits = []
            for p_ in range(0, k + 1):
                for i in range(0, min(p_, k - 1) + 1):
                    y = n if (i == p_ or (p_ == k and i == k - 1)) else (n, base + 3 * i)
                    orbits.append([3, y] if p_ == k else [3, y, base + 3 * p_])
            if k == 4 and n == 14:
                print("   k=4 [3,3k+2]: beyond truncation K=48, not computed"); continue
            inv = [a_b(c)[:2] + (tuple(sorted(local_algebra([series(x) for x in c]).rows))[:12],) for c in orbits]
            distinct = len(set(inv)) == len(inv)
            print(f"   k={k} {name}: {len(orbits):2} orbits (k(k+3)/2 = {k * (k + 3) // 2}); "
                  f"pairwise separated by (a, b, semigroup): {distinct}")

  if "C" in parts:
    print("\nC. Arnold Thm 4 with i > p: (3,(10,14),11) reduces to (3,10,11) since t^14 = x*z on the curve.")
    print("   computed codim of (3,(10,14),11):", "%dN-%d" % a_b([3, (10, 14), 11])[:2],
          "; Arnold's formula at k=3,p=0,i=1 gives 6N-3; E_{18,0,0}=(3,10,11) has", "%dN-%d" % a_b([3, 10, 11])[:2])
