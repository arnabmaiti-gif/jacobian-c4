"""Degree-4 analogue of the cubic-factorization Jacobian counterexample.

L(t) = x t + beta (linear), Q(t) = g3 t^3 + g2 t^2 + g1 t + g0 (cubic).
Product quartic c4..c0, resultant rho = Res(L,Q) = x^3 Q(-beta/x)
                                     = -beta^3 g3 + x beta^2 g2 - x^2 beta g1 + x^3 g0.
Slice: c3 = 1, rho = -1.  Chart branch beta = 1 + x*y1.
"""
import sympy as sp

x, y1, y2, y3 = sp.symbols("x y1 y2 y3")
X, B, G3, G2, G1, G0 = sp.symbols("X B G3 G2 G1 G0")

# ---------- Step 1: coefficient-resultant Jacobian identity (generic coords) ----------
c4g, c3g, c2g, c1g, c0g = X * G3, X * G2 + B * G3, X * G1 + B * G2, X * G0 + B * G1, B * G0
rho_g = -B**3 * G3 + X * B**2 * G2 - X**2 * B * G1 + X**3 * G0
# sanity: rho_g equals sympy's resultant up to sign convention
L = X * sp.Symbol("t") + B
Q = G3 * sp.Symbol("t") ** 3 + G2 * sp.Symbol("t") ** 2 + G1 * sp.Symbol("t") + G0
res_builtin = sp.resultant(L, Q, sp.Symbol("t"))
print("rho matches sympy resultant:", sp.expand(rho_g - res_builtin) == 0,
      "| up to sign:", sp.expand(rho_g + res_builtin) == 0)

M = sp.Matrix([c4g, c3g, c2g, c1g, c0g, rho_g]).jacobian([X, B, G3, G2, G1, G0])
detM = sp.factor(M.det())
print("det d(c4,c3,c2,c1,c0,rho)/d(6 params) =", detM)

# ---------- Step 2: polynomial chart of the slice {c3 = 1, rho = -1} ----------
beta = 1 + x * y1
g = -2 * y1 + x * y2                      # gamma3 = 1 + x*g ; the -2*y1 cancels the x^1 obstruction
gamma3 = sp.expand(1 + x * g)
num = sp.expand(1 - beta * gamma3)        # c3 = x*gamma2 + beta*gamma3 = 1  =>  gamma2 = (1-beta*gamma3)/x
q2, r2 = sp.div(num, x, x)
assert r2 == 0, ("gamma2 not polynomial", r2)
gamma2 = sp.expand(q2)

w = sp.expand(2 * beta**2 * gamma2 - y1 * (beta + 1))   # rho=-1 <=> x*beta*gamma1 - x^2*gamma0 = w
w1, rw = sp.div(w, x, x)
assert rw == 0, ("w not divisible by x", rw)
gamma1 = sp.expand(w1 + x * y3)
gamma0 = sp.expand(y1 * w1 + beta * y3)

# hard checks: slice equations hold identically
c3 = sp.expand(x * gamma2 + beta * gamma3)
rho = sp.expand(-beta**3 * gamma3 + x * beta**2 * gamma2 - x**2 * beta * gamma1 + x**3 * gamma0)
print("c3 == 1:", c3 == 1)
print("rho == -1:", rho == -1)

# ---------- Step 3: the map = retained coefficients ----------
c4 = sp.expand(X * G3).subs({X: x, G3: gamma3})
c2 = sp.expand(x * gamma1 + beta * gamma2)
c1 = sp.expand(x * gamma0 + beta * gamma1)
c0 = sp.expand(beta * gamma0)
F = sp.Matrix([c4, c2, c1, c0])
detJ = sp.expand(F.jacobian([x, y1, y2, y3]).det())
print("det JF =", detJ)

for name, comp in zip("F1 F2 F3 F4".split(), F):
    print(name, "=", comp)
