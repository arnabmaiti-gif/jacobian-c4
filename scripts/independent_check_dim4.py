"""Independent check of the C^4 counterexample, from the printed polynomials only."""
import sympy as sp

x, y1, y2, y3 = sp.symbols("x y1 y2 y3")

F = sp.Matrix([
    x**3*y2 - 2*x**2*y1 + x,
    -2*x**4*y1**3*y2 + 4*x**3*y1**4 - 7*x**3*y1**2*y2 + 12*x**2*y1**3 - 8*x**2*y1*y2
        + x**2*y3 + 10*x*y1**2 - 3*x*y2 + y1,
    -4*x**4*y1**4*y2 + 8*x**3*y1**5 - 14*x**3*y1**3*y2 + 24*x**2*y1**4 - 18*x**2*y1**2*y2
        + 2*x**2*y1*y3 + 24*x*y1**3 - 10*x*y1*y2 + 2*x*y3 + 7*y1**2 - 2*y2,
    -2*x**4*y1**5*y2 + 4*x**3*y1**6 - 8*x**3*y1**4*y2 + 14*x**2*y1**5 - 12*x**2*y1**3*y2
        + x**2*y1**2*y3 + 17*x*y1**4 - 8*x*y1**2*y2 + 2*x*y1*y3 + 7*y1**3 - 2*y1*y2 + y3,
])

det = sp.expand(F.jacobian([x, y1, y2, y3]).det())
print("det JF =", det)
assert det == -2

pts = [(1, -1, -3, -1), (-1, 1, -3, 1), (0, 0, sp.Rational(1, 2), 0)]
vals = [F.subs(dict(zip((x, y1, y2, y3), p))) for p in pts]
for p, v in zip(pts, vals):
    print(f"F{p} =", v.T)
assert vals[0] == vals[1] == vals[2]
assert len(set(pts)) == 3
print("VERIFIED: constant Jacobian -2, three distinct rational points collide.")
