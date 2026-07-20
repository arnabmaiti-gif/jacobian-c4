"""One-file exact verification of the degree-4 Jacobian-conjecture counterexample family member.

Construction: L(t) = x t + beta,  Q(t) = g3 t^3 + g2 t^2 + g1 t + g0,
slice {c3 = 1, rho = -1} of the coefficient+resultant map, polynomial chart branch
beta = 1 + x y1,  g3 = 1 - 2 x y1 + x^2 y2,  g1 = w1 + x y3  (w1 forced),  g0 = y1 w1 + beta y3.
Output map F = (c4, c2, c1, c0): A^4 -> A^4.
"""
import sympy as sp

x, y1, y2, y3, t = sp.symbols("x y1 y2 y3 t")

# ---- chart ----
beta = 1 + x * y1
gamma3 = sp.expand(1 + x * (-2 * y1 + x * y2))
gamma2, r = sp.div(sp.expand(1 - beta * gamma3), x, x)
assert r == 0
w1, r = sp.div(sp.expand(2 * beta**2 * gamma2 - y1 * (beta + 1)), x, x)
assert r == 0
gamma1 = sp.expand(w1 + x * y3)
gamma0 = sp.expand(y1 * w1 + beta * y3)

# ---- slice identities ----
assert sp.expand(x * gamma2 + beta * gamma3) == 1                       # c3 == 1
rho = -beta**3 * gamma3 + x * beta**2 * gamma2 - x**2 * beta * gamma1 + x**3 * gamma0
assert sp.expand(rho) == -1                                             # resultant == -1
print("slice identities: c3 = 1 and Res(L,Q) = -1 hold identically")

# ---- the map ----
F = sp.Matrix([
    sp.expand(x * gamma3),                    # c4
    sp.expand(x * gamma1 + beta * gamma2),    # c2
    sp.expand(x * gamma0 + beta * gamma1),    # c1
    sp.expand(beta * gamma0),                 # c0
])
print("\nF components (A^4 -> A^4):")
for i, comp in enumerate(F, 1):
    print(f"  F{i} =", comp)

detJ = sp.expand(F.jacobian([x, y1, y2, y3]).det())
print("\ndet JF =", detJ)
assert detJ == -2

# ---- the 7-point fiber over the quartic t^3 - t : target (0,0,-1,0) ----
s2 = sp.sqrt(2) * sp.I
fiber = [
    (1, -1, -3, -1),
    (-1, 1, -3, 1),
    (0, 0, sp.Rational(1, 2), 0),
    (s2 / 2, -1 + s2, 6 + 2 * s2, -2 - s2),
    (-s2 / 2, -1 - s2, 6 - 2 * s2, -2 + s2),
    (s2 / 2, 1 + s2, 6 - 2 * s2, 2 - s2),
    (-s2 / 2, 1 - s2, 6 + 2 * s2, 2 + s2),
]
target = sp.Matrix([0, 0, -1, 0])
print("\nfiber over (c4,c2,c1,c0) = (0,0,-1,0):")
for p in fiber:
    val = sp.expand(F.subs(dict(zip((x, y1, y2, y3), p))))
    assert sp.simplify(val - target) == sp.zeros(4, 1), p
    print("  F", tuple(p), "= (0,0,-1,0)  [exact]")
assert len(set(map(tuple, fiber))) == 7
print("7 distinct points, exact arithmetic. F is 7:1 over this target => not injective.")

# ---- general-n pattern: coefficient-resultant Jacobian det for L(1) x Q(m) ----
print("\ncoefficient-resultant Jacobian identity, L linear x Q of degree m:")
X, B = sp.symbols("X B")
for m in (2, 3, 4):
    gs = sp.symbols(f"g0:{m + 1}")            # Q = gm t^m + ... + g0
    Q = sum(gs[i] * t**i for i in range(m + 1))
    Lp = X * t + B
    prod = sp.expand(Lp * Q)
    cs = [prod.coeff(t, k) for k in range(m + 1, -1, -1)]
    rho_m = sp.expand((-1)**m * X**m * Q.subs(t, -B / X) * sp.Rational(1))
    rho_m = sp.expand(sp.cancel(X**m * Q.subs(t, -B / X)))
    Mj = sp.Matrix(cs + [rho_m]).jacobian([X, B] + list(gs))
    d = sp.factor(Mj.det())
    print(f"  m={m} (degree n={m + 1}):  det = {d}")
