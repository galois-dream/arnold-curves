"""Running times of the note's three algorithms on growing inputs (are they practical, not just stated?).

Written 2026-10-15 23:53 CEST.  [Ours]  One file because it times all three algorithms side by side; the
algorithms themselves live in codim_check.codimension (Alg. 1), normalize.normalize (Alg. 2) and
graded_coker.eliminate (Alg. 3).  Pure Python, exact rational arithmetic, one core; prints seconds.
Run: .venv/bin/python computations/timing.py          (the examples quoted in the note)
     .venv/bin/python computations/timing.py sweep    (growth rates: time and peak memory against size, with
                                                       fitted log-log slopes; added 2026-10-16 00:09 CEST)
"""
import math
import sys
import time
import tracemalloc

import codim_check
import normalize as nz
import graded_coker
from codim_check import codimension
from graded_coker import curve_moves, eliminate
from normalize import decide, normalize, t


def clock(f, *args, **kw):
    t0 = time.perf_counter()
    out = f(*args, **kw)
    return out, time.perf_counter() - t0


def measure(f, *args):
    """(seconds, peak MiB): time from a plain run, memory from a second run under tracemalloc."""
    _, s = clock(f, *args)
    tracemalloc.start()
    f(*args)
    peak = tracemalloc.get_traced_memory()[1] / 2**20
    tracemalloc.stop()
    return s, peak


def slope(xs, ys):
    """Least-squares slope of log y against log x, over the larger half of the sizes."""
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if y > 0][len(xs) // 2 - 1:]
    mx, my = sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
    return sum((x - mx) * (y - my) for x, y in pts) / sum((x - mx) ** 2 for x, _ in pts)


def sweep():
    def report(name, size_name, rows, predicted):
        xs = [r[0] for r in rows]
        for r in rows:
            print(f"    {size_name}={r[0]:<5} {r[1]:9.3f} s  {r[2]:9.2f} MiB" + (f"   ({r[3]})" if len(r) > 3 else ""))
        print(f"  {name}: time ~ {size_name}^{slope(xs, [r[1] for r in rows]):.1f}, "
              f"memory ~ {size_name}^{slope(xs, [r[2] for r in rows]):.1f}   ({predicted})")

    for label, curve in (("(t^3, t^n), sparse", lambda n: [3, n]),
                         ("(t^3+t^4+t^5, t^n+t^(n+1)+t^(n+2)), dense", lambda n: [(3, 4, 5), (n, n + 1, n + 2)])):
        rows = []
        for n in (5, 10, 14, 20, 28, 40):
            s, mb = measure(codimension, curve(n))
            rows.append((n, s, mb, f"K={codimension(curve(n))[2]}"))
        report(f"Alg. 1 on {label}", "n", rows, "worst case time n^3, memory n^2 (K grows like n)")

    rows = []
    for n in (4, 7, 10, 14, 20, 28, 40):
        nz.ORD = 3 * n + 6
        codim_check.K = max(48, 4 * n + 8)  # normalize() reads the conductor through a_b
        s, mb = measure(nz.normalize, 3, n, t**(2 * n - 3))
        rows.append((n, s, mb))
    nz.ORD, codim_check.K = 60, 48
    report("Alg. 2 on (t^3, t^n + t^(2n-3)), cut off at 3n+6", "n", rows, "worst case time n^2, memory n^2: each substitution is expanded before it is cut off")

    for curve, lead, pred in (([4, (6, 7)], [4, 6], "worst case time K^4, memory K^3"),
                              ([4, (6, 7), 9], [4, 6, 9], "worst case time K^5, memory K^4")):
        rows = []
        for K in (32, 48, 64, 96, 128):
            codim_check.K = graded_coker.K = K
            cols, degs, rws = curve_moves(curve, lead)
            s, mb = measure(lambda: eliminate(*curve_moves(curve, lead)))
            rows.append((K, s, mb, f"{len(cols)} x {len(rws)}"))
        codim_check.K = graded_coker.K = 48
        report(f"Alg. 3 on {curve} (building the matrix + elimination)", "K", rows, pred)
    codim_check.K = graded_coker.K = 96
    cols, degs, rws = curve_moves([4, (6, 7), 9], [4, 6, 9])
    piv = {}
    eliminate(cols, degs, rws, piv)
    bits = lambda xs: max(max(abs(x.numerator).bit_length(), x.denominator.bit_length()) for x in xs)
    print(f"  coefficient size, Alg. 3 on [4, (6, 7), 9] at K=96: inputs <= {bits([graded_coker.Fraction(c) for col in cols for c in col.values()])} bits,"
          f" pivots <= {bits([x for v in piv.values() for x in v.values()])} bits")
    codim_check.K = graded_coker.K = 48


if __name__ == "__main__":
    if sys.argv[1:] == ["sweep"]:
        sweep()
        sys.exit()
    print("Algorithm 1 (codimension aN - b, K chosen automatically)")
    for curve in ([2, 3], [3, 4], [3, 31], [4, (6, 7)], [4, (6, 21)], [4, 7, (9, 10)], [3, (13, 23), 17],
                  [5, 6, 7, 8, 9], [6, 7, 8, 9, 10, 11], [5, (6, 8), 7], [7, 11]):
        (a, b, K), s = clock(codimension, curve)
        print(f"  {str(curve):24} {a}N-{b:<4} K={K:<4} {s:7.2f} s")

    print("Algorithm 2 (term removal) on (t^3, t^n + t^(2n-3)), then the full Algorithm 2 (invariant test first)")
    for n in (4, 5, 7, 8, 10, 13, 16):
        (log, rest, cond), s = clock(normalize, 3, n, t**(2 * n - 3))
        print(f"  n={n:<3} {len(log)} steps, stopped at order {rest} (conductor {cond})   {s:6.2f} s")
    for m, n, extra in ((3, 7, t**8), (3, 7, t**11), (3, 7, t**8 + t**11), (4, 9, t**10), (4, 9, t**23)):
        (verdict, why), s = clock(decide, m, n, extra)
        print(f"  decide (t^{m}, t^{n} + {extra}): {verdict} ({why})   {s:6.2f} s")

    print("Algorithm 3 (ordered elimination) on the curve's moves, as the truncation order K grows")
    for curve, lead in (([4, (6, 7)], [4, 6]), ([4, (6, 7), 9], [4, 6, 9])):
        for K in (32, 48, 64, 96):
            codim_check.K = graded_coker.K = K
            cols, degs, rows = curve_moves(curve, lead)
            (_, free), s = clock(eliminate, cols, degs, rows)
            print(f"  {str(curve):16} K={K:<3} {len(cols):5} columns x {len(rows):4} rows   {s:6.2f} s")
