# arnold-curves

Supporting code for the note *Gaps, Curves, and Arnold's Spectral Sequence* (draft, October 2026).
Every number quoted in the note is printed by one of these scripts. All linear algebra uses exact rational
arithmetic; symbolic substitutions use `sympy`.

Written 2026-10-15 15:25 CEST.

| script | what it reproduces |
|---|---|
| `computations/codim_check.py` | Algorithm 1: Arnold's codimension formula aN − b (gaps of the semigroup, rank of one matrix), with the cutoff chosen automatically; checked against the codimensions printed in Arnold (1999), including the e_k / f_k bracket check in his Theorem 5 |
| `computations/note_examples.py` | the small worked examples: (t^2, t^3+t^4), (t^3, t^4+t^5), the count b for (t^3, t^4), Sylvester duality for 3 and 5, the rolling wheel, the codimensions of (2,3), (3,4,5), (3,4) |
| `computations/normalize.py` | Algorithm 2: explicit substitution proofs that the Frobenius ("last gap") term can be removed, e.g. (t^3, t^n + t^(2n−3)) ~ (t^3, t^n); `decide()` adds the codimension test and a step budget |
| `computations/cm_3_10_orbits.py` | the nine curve classes with symbol [3,10] and their codimensions against Arnold's E_18 formula |
| `computations/moduli_test.py` | codimension of a symbol class vs a generic orbit (where moduli appear); the count k(k+3)/2; Arnold's Theorem 4 index check |
| `computations/graded_coker.py` | the graded cokernel and the first differentials (equation, Euler and Hamiltonian fields) for (t^4, t^6 + t^7); Algorithm 3 (ordered elimination), checked against the graded cokernel |
| `computations/timing.py` | running times of Algorithms 1–3 on the note's examples; `timing.py sweep` measures time and peak memory against input size (about 30 s) |
| `computations/gibson_hobbs_check.py` | Gibson–Hobbs (1993) Tables 1–2 against aN − b at N = 3 (92 normal forms), and the entries their invariants cannot separate (the Frobenius pairs) |
| `computations/plane_branches.py` | plane curves: the two-move removal of the Frobenius term for all coprime m < n (Lemma 2); Zariski's theorem by elimination and the Euler move on 158 curves (Proposition 3); Arnold's b against Hefez–Hernandes's values of differentials (b = a − 1 + #(Λ∖Γ)) and the codimension at N = 3 against the Tjurina number (about 15 s) |
| `computations/codewords.py` | code-word counts for the Semple tower strata, the RVT ↔ code-word bijection, closed-form generating functions |

## Running

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python computations/codim_check.py      # each script runs on its own
.venv/bin/python computations/moduli_test.py A    # part A takes about two minutes
```

Python 3.10 or later. Each script's docstring states what it checks and labels its sources:
[Arnold] for Arnold's paper (page references are to the Russian edition, Trudy MIAN 226, 1999),
[Lit] for other literature, and [Ours] for new reconstructions.

## License

MIT; see `LICENSE`.
