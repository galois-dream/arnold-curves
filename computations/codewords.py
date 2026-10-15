"""Reproduce the code-word combinatorics of Castro–Colley–Kennedy–Shanbrom, arXiv:1606.07931 (CCKS), and
cross-check them against the RVT spelling rules of Castro–Howard–Shanbrom, arXiv:1512.00107 (CHS).

Written 2026-10-15 15:18 CEST (corrected from a guessed 15:33).  [Ours]
CCKS §3 rules (base dimension m): symbols V_A, A a finite set of integers > 1, R = V_{} ;
  (1) first symbol R; (2) after V_A at position j-1, the symbol V_B at position j has B ⊆ A ∪ {j};  (3) |A| < m.
  [We read rule (2) "B ⊆ A, or B ⊆ A ∪ {j}" as B ⊆ A ∪ {j}.]
CCKS §7: N(k,r) recursion; N(k,r) = c(k,r+1) (unsigned Stirling, first kind) when k <= m.
CCKS Thm 4.2: codim of the stratum = sum of |A_j|.
CHS Thm 2 (= CCKS §3 RVT rules, m = 3): R -> R|V ; V,T1 -> R|V|T1|L1 ; T2 -> R|V|T2|L3 ; L1,L2,L3 -> anything.
CCKS translation tables (p.6) map RVT words to code words; we check the map is a bijection that preserves codim
(V,T1,T2 -> 1 ; L1,L2,L3 -> 2, as in CM2012 where L-classes have codim one more than V-classes).
Run: .venv/bin/python computations/codewords.py
"""
from itertools import combinations
from math import factorial


def code_words(k, m):
    words = [((),)]  # tuple of frozensets; position 1 is R = empty set
    words = [(frozenset(),)]
    for j in range(2, k + 1):
        new = []
        for w in words:
            pool = sorted(w[-1] | {j})
            for r in range(0, min(len(pool), m - 1) + 1):
                for B in combinations(pool, r):
                    new.append(w + (frozenset(B),))
        words = new
    return words


def stirling1(n, k):
    if n == k:
        return 1
    if n == 0 or k == 0:
        return 0
    return (n - 1) * stirling1(n - 1, k) + stirling1(n - 1, k - 1)


NEXT = {"R": "RV", "V": ["R", "V", "T1", "L1"], "T1": ["R", "V", "T1", "L1"], "T2": ["R", "V", "T2", "L3"],
        "L1": ["R", "V", "T1", "T2", "L1", "L2", "L3"], "L2": ["R", "V", "T1", "T2", "L1", "L2", "L3"],
        "L3": ["R", "V", "T1", "T2", "L1", "L2", "L3"]}


def rvt_words(k):
    words = [["R"]]
    for _ in range(k - 1):
        words = [w + [x] for w in words for x in (["R", "V"] if w[-1] == "R" else NEXT[w[-1]])]
    return words


def translate(w):
    """CCKS p.6 tables: RVT word -> code word (list of sets). The table used depends on the previous code symbol."""
    out = [frozenset()]
    for k in range(2, len(w) + 1):
        prev, sym = out[-1], w[k - 1]
        if not prev:  # after R
            tab = {"R": set(), "V": {k}}
        elif len(prev) == 1:
            (j,) = tuple(prev)
            after_T2 = w[k - 2] == "T2"
            tab = ({"R": set(), "V": {k}, "T2": {j}, "L3": {j, k}} if after_T2 else
                   {"R": set(), "V": {k}, "T1": {j}, "L1": {j, k}})
        else:
            i, j = sorted(prev)
            tab = {"R": set(), "V": {k}, "T1": {j}, "T2": {i}, "L1": {j, k}, "L2": {i, j}, "L3": {i, k}}
        if sym not in tab:
            return None
        out.append(frozenset(tab[sym]))
    return tuple(out)


if __name__ == "__main__":
    print("CCKS §7: number of code words of length k (base dimension m), and Stirling check for k <= m")
    for m in (2, 3, 4):
        counts = [len(code_words(k, m)) for k in range(1, 8)]
        stir_ok = all(sum(1 for w in code_words(k, m) if len(w[-1]) == r) == stirling1(k, r + 1)
                      for k in range(1, m + 1) for r in range(0, k))
        print(f"  m={m}: lengths 1..7 -> {counts};  N(k,r)=c(k,r+1) for k<=m: {stir_ok};  k! = "
              f"{[factorial(k) for k in range(1, m + 1)]}")
    print("\nCross-check with CHS RVT spelling rules (m = 3):")
    for k in range(1, 8):
        rv = rvt_words(k)
        cw = set(code_words(k, 3))
        img = [translate(w) for w in rv]
        bad = [w for w, t in zip(rv, img) if t is None or t not in cw]
        bij = len(set(i for i in img if i)) == len(rv) == len(cw) and not bad
        codim_ok = all(sum(1 if s in ("V", "T1", "T2") else 2 if s.startswith("L") else 0 for s in w) ==
                       sum(len(A) for A in t) for w, t in zip(rv, img) if t)
        print(f"  k={k}: RVT words {len(rv):4}, code words {len(cw):4}, translation is a bijection: {bij}, "
              f"codim preserved: {codim_ok}" + (f"   e.g. untranslatable {''.join(bad[0])}" if bad else ""))
    ex = translate(["R", "V", "T1", "V", "V", "T1", "T1", "T1", "L1"])
    print("\nCCKS p.7 example RVT1VVT1T1T1L1 ->", " ".join("R" if not A else "V" + "".join(map(str, sorted(A))) for A in ex),
          "(paper: R V2 V2 V4 V5 V5 V5 V5 V59)")

    # [Ours] closed forms: CCKS §7 recursion as a transfer matrix M[r][i] = C(i+1, r), r, i in 0..m-1
    import sympy as sp
    print("\nClosed forms from the CCKS §7 transfer matrix M[r][i] = binom(i+1, r):")
    x = sp.symbols("x")
    for m in (2, 3, 4):
        M = sp.Matrix(m, m, lambda r, i: sp.binomial(i + 1, r))
        e0 = sp.Matrix([1] + [0] * (m - 1))
        ones = sp.ones(1, m)
        seq = [(ones * M**(k - 1) * e0)[0] for k in range(1, 9)]
        enum = [len(code_words(k, m)) for k in range(1, 8)]
        gf = sp.simplify((ones * (sp.eye(m) - x * M).inv() * e0)[0])
        cp = M.charpoly(x).as_expr()
        lam = max(sp.Poly(cp, x).nroots(), key=lambda z: abs(z))
        print(f"  m={m}: counts {seq[:7]} (enumeration agrees: {seq[:7] == enum});  sum_k count_k x^(k-1) = {sp.factor(gf)};"
              f"  growth ~ {sp.N(lam, 6)}^k")
