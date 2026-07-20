"""Direct construction + exact verification of a 7-point fiber of the degree-4 map.

Fiber over the degenerate quartic  t^3 - t  <->  target (c4,c2,c1,c0) = (0,0,-1,0).
Candidates: 3 affine roots x 2 resultant scalings + 1 infinity sheet on the chart branch.
"""
import sympy as sp

x, y1, y2, y3 = sp.symbols("x y1 y2 y3")

# chart (as constructed): beta, gamma3, gamma2, gamma1, gamma0 as polys in (x,y1,y2,y3)
beta = 1 + x * y1
gamma3 = sp.expand(1 + x * (-2 * y1 + x * y2))
gamma2 = sp.expand(sp.div(sp.expand(1 - beta * gamma3), x, x)[0])
w = sp.expand(2 * beta**2 * gamma2 - y1 * (beta + 1))
w1 = sp.div(w, x, x)[0]
gamma1 = sp.expand(w1 + x * y3)
gamma0 = sp.expand(y1 * w1 + beta * y3)

F = sp.Matrix([
    sp.expand(x * gamma3),
    sp.expand(x * gamma1 + beta * gamma2),
    sp.expand(x * gamma0 + beta * gamma1),
    sp.expand(beta * gamma0),
])
assert sp.expand(F.jacobian([x, y1, y2, y3]).det()) == -2

lam = sp.I / sp.sqrt(2)   # lambda^2 = -1/2, the scaling putting the theta=+-1 sheets on rho=-1

# candidate Y-points as (x, beta, g3, g2, g1, g0)
cands = {
    "A (L=t)":        (1, 0, 0, 1, 0, -1),
    "B (L=-t)":       (-1, 0, 0, -1, 0, 1),
    "C (infinity)":   (0, 1, 1, 0, -1, 0),
    "D+ (L~t-1)":     (lam, -lam, 0, 1/lam, 1/lam, 0),
    "D- (L~t-1)":     (-lam, lam, 0, -1/lam, -1/lam, 0),
    "E+ (L~t+1)":     (lam, lam, 0, 1/lam, -1/lam, 0),
    "E- (L~t+1)":     (-lam, -lam, 0, -1/lam, 1/lam, 0),
}

target = sp.Matrix([0, 0, -1, 0])
pts = []
for name, cand in cands.items():
    xv, bv, g3v, g2v, g1v, g0v = [sp.sympify(v) for v in cand]
    if xv != 0:
        Y1 = sp.radsimp((bv - 1) / xv)
        Y2 = sp.radsimp((g3v - 1 + 2 * xv * Y1) / xv**2)
        Y3 = sp.radsimp((g1v - w1.subs({x: xv, y1: Y1, y2: Y2})) / xv)
    else:  # infinity sheet: x=0, beta=1 forces y1 via gamma2, y2 via gamma1, y3 via gamma0
        Y1 = sp.solve(gamma2.subs(x, 0) - g2v, y1)[0]
        Y2 = sp.solve(gamma1.subs({x: 0, y1: Y1}) - g1v, y2)[0]
        Y3 = sp.solve(gamma0.subs({x: 0, y1: Y1, y2: Y2}) - g0v, y3)[0]
    p = {x: xv, y1: Y1, y2: Y2, y3: Y3}
    val = sp.simplify(F.subs(p))
    ok = sp.simplify(val - target) == sp.zeros(4, 1)
    pts.append((sp.nsimplify(xv), Y1, Y2, Y3))
    print(f"{name}: point = ({xv}, {Y1}, {Y2}, {Y3})  F(point)==target: {ok}")

distinct = len({tuple(sp.simplify(c) for c in p) for p in pts})
print("distinct points:", distinct)
