"""Greedy term-killing normalizer for plane curves (t^m, t^n + ...): proves A-equivalences by explicit substitutions.

Written 2026-10-15 ~15:05 CEST (approx.; committed 15:16).  [Ours], following the scheme of [Arnold] p.34:
  repeat: take the lowest-order extra term c t^e (over both slots);
    if e is in the semigroup of the target and e > leading exponent of that slot -> left change x_i -= c P(x, y),
       where P is a monomial x^alpha y^beta with alpha*m + beta*n = e (other than x_i itself);
    else (e a gap) -> right change t -> t (1 - c/lead_i t^(e - lead_i)), which kills it at first order;
  stop when every remaining extra term has order >= the conductor (then left changes remove them all).
Every step is an explicit invertible substitution, carried out exactly in sympy (truncated at order ORD).
Run: .venv/bin/python computations/normalize.py
"""
import sympy as sp

from codim_check import a_b

t = sp.symbols("t")
ORD = 60


def semigroup(m, n, upto):
    return {a * m + b * n for a in range(upto // m + 1) for b in range(upto // n + 1) if 0 < a * m + b * n <= upto}


def trunc(e):
    e = sp.expand(e)
    return sum(e.coeff(t, k) * t**k for k in range(ORD))


def normalize(m, n, extra, max_steps=200):
    """Start from (t^m, t^n + extra); return (steps, remaining_min_order, conductor)."""
    cond = a_b([m, n])[2]
    G = semigroup(m, n, ORD)
    lead = [m, n]
    f = [t**m, trunc(t**n + extra)]
    log = []
    for _ in range(max_steps):
        extras = []
        for i in range(2):
            r = sp.Poly(sp.expand(f[i] - t**lead[i]), t)
            extras += [(mon[0], i, c) for mon, c in r.terms()] if r.as_expr() != 0 else []
        if not extras:
            return log, None, cond
        e, i, c = min(extras)
        if e >= cond:
            return log, e, cond
        if e in G and e > lead[i]:
            # find P = x^a y^b with a*m + b*n = e, not x_i itself
            a, b = next((a, b) for a in range(e // m + 1) for b in range(e // n + 1)
                        if a * m + b * n == e and (a, b) != ((1, 0) if i == 0 else (0, 1)))
            P = f[0]**a * f[1]**b
            lc = sp.Poly(sp.expand(P), t).coeff_monomial(t**e)
            f[i] = trunc(f[i] - (c / lc) * P)
            log.append(f"left  x{i + 1} -= ({c / lc})*x^{a} y^{b}   (kills t^{e})")
        else:
            r = e - lead[i]
            s = -c / lead[i]
            f = [trunc(fi.subs(t, t * (1 + s * t**r))) for fi in f]
            log.append(f"right t -> t(1 + ({s}) t^{r})          (kills t^{e} in slot {i + 1})")
    raise RuntimeError("did not terminate")



def decide(m, n, extra, max_steps=12):
    """Algorithm 2 of the note as stated there (added 2026-10-15 23:53 CEST, corrected from a guessed 23:58): first the invariant test, then
    the term-removal loop with a step budget.  Returns ('not equivalent' | 'equivalent' | 'undecided', detail)."""
    from codim_check import codimension
    terms = sp.Poly(sp.expand(t**n + extra), t).terms()
    curve = [m, {e[0]: c for e, c in terms}]
    try:
        lead = codimension([m, n])[:2]
    except ValueError:
        lead = None  # leading curve infinitely degenerate
    if codimension(curve)[:2] != lead:
        return "not equivalent", f"codimensions differ: {codimension(curve)[:2]} vs {lead}"
    try:
        log, rest, cond = normalize(m, n, extra, max_steps=max_steps)
        return "equivalent", f"{len(log)} explicit steps"
    except RuntimeError:
        return "undecided", f"no answer within {max_steps} steps"

CASES = [  # (m, n, extra, citation)
    (3, 4, t**5, "Arnold E_{6,1,0} boundary: (3,(4,5)) ~ (3,4)"),
    (3, 5, t**7, "Arnold E_{8,1,0}: (3,(5,7)) ~ (3,5); Castro-Howard 1107.4145 p.17 'messy calculation'"),
    (3, 5, -t**7, "same, minus sign"),
    (3, 7, t**11, "Arnold E_{12,2,1} = (3,7): (3,(7,11)) ~ (3,7); CM2012 Remark 3.33"),
    (3, 8, t**13, "Arnold E_{14,2,1} = (3,8): (3,(8,13)) ~ (3,8)"),
    (3, 10, t**17, "Arnold E_{18,3,2} = (3,10): (3,(10,17)) ~ (3,10); CM arXiv v1 p.12"),
    (3, 13, t**23, "Arnold E_{24,4,3} = (3,13): (3,(13,23)) ~ (3,13)"),
    # Frobenius term F = (m-1)(n-1)-1 of <m,n> for m >= 4 (added 2026-10-15, committed 15:16 CEST)
    (4, 5, t**11, "Frobenius term of <4,5>"),
    (4, 7, t**17, "Frobenius term of <4,7>"),
    (4, 9, t**23, "Frobenius term of <4,9>"),
    (5, 6, t**19, "Frobenius term of <5,6>"),
    (5, 7, t**23, "Frobenius term of <5,7>"),
]

if __name__ == "__main__":
    for m, n, extra, note in CASES:
        log, rest, cond = normalize(m, n, extra)
        ok = rest is None or rest >= cond
        print(f"(t^{m}, t^{n} + {extra}) ~ (t^{m}, t^{n}): {'PROVED' if ok else 'NOT SHOWN'} "
              f"in {len(log)} explicit steps; remaining terms of order >= {rest} (conductor {cond}).  [{note}]")
        for line in log:
            print("     ", line)
