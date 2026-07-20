# An essentially four-dimensional counterexample to the Jacobian conjecture

**Webpage: https://arnabmaiti-gif.github.io/jacobian-c4/**

An explicit degree-10 polynomial map F: C⁴ → C⁴ with integer coefficients,
constant Jacobian determinant −2, and a verified seven-point fiber (including a
three-point rational collision). Its generic fiber has ≥ 7 points, so it is not a
re-coordinatization or identity-padding of Alpöge's map, and the construction
provably cannot be carried out in C³ — to our knowledge the first counterexample
*constructed natively* in dimension 4. (Whether it is equivalent to a padding of a
same-fiber-degree member of the jacobianfun C³ family is not settled by fiber
counts; see §4 of the webpage.)

Built by running the cubic-factor/resultant construction behind Alpöge's
July 19, 2026 counterexample in C³ one degree up (quartic = linear × cubic,
slice {c₃ = 1, Res = −1}).

## Verify

All checks are exact rational arithmetic (SymPy), no floating point:

```bash
cd scripts
uv run --with sympy python independent_check_dim4.py  # theorem only, from the displayed polynomials
uv run --with sympy python verify_all.py              # full construction + all 7 fiber points
```

- `independent_check_dim4.py` — minimal certificate: det JF ≡ −2 and the rational 3-collision.
- `construct.py` — chart construction and the coefficient–resultant Jacobian identity.
- `fiber_direct.py` — derives the 7 fiber points from the factorization classification.
- `verify_all.py` — consolidated assert-everything run.

## Credits

Original C³ counterexample: Levent Alpöge with Claude Fable 5 (question posed by
Akhil Mathew). Derivation this generalizes:
[aaronlou.com/jacobian_counterexample_derivation.pdf](https://aaronlou.com/jacobian_counterexample_derivation.pdf).
C³ family of every generic fiber degree: [jacobianfun.org](https://jacobianfun.org).
This note: Arnab Maiti and Esha Manideep, with the construction and verification
carried out with Claude (Fable 5) agents. Not peer-reviewed; corrections welcome via issues.
