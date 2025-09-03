import sympy as sp
import argparse
from math import ceil

# Parse arguments
parser = argparse.ArgumentParser(description="Generate inlet polynomial velocity profile with given flow rate")
parser.add_argument("--polynomial-degree", type=int, default=2, help="Degree k of the polynomial shape (k>=0)")
parser.add_argument("--shape", type=str, default="rectangular", choices=["rectangular", "circular"],
                    help="Inlet cross-section shape")
parser.add_argument("--flow-rate", type=float, default=1.0, help="Desired volumetric flow rate Q")
parser.add_argument("--center", type=float, nargs="+", default=[0.0, 0.0],
                    help="Center coordinates [x0, y0]")
parser.add_argument("--length", type=float, nargs="+", default=[1.0, 1.0],
                    help="Rectangular: lengths [ax, ay]. Circular: diameter [D]")
parser.add_argument("--show", action="store_true", help="Show a 1D centerline plot (rectangular only)")
args = parser.parse_args()

# Symbols
x, y = sp.symbols('x y', real=True)
x0, y0 = map(float, args.center)
k = int(args.polynomial_degree)
Q = float(args.flow_rate)

if k < 0:
    raise ValueError("polynomial-degree must be >= 0")

# Odd degrees lead to negative velocity
if k % 2 == 1:
    raise ValueError("Odd polynomial degrees are not supported.")

# Build polynomial shape p(x,y) and integrate it over the inlet
if args.shape == "rectangular":
    # length = [ax, ay] are lengths; domain is [x0-ax, x0+ax] x [y0-ay, y0+ay]
    if len(args.length) != 2:
        raise ValueError("For rectangular, --length must be two numbers: [ax ay] (lengths).")
    half_lengths = [l / 2 for l in args.length]
    ax, ay = half_lengths
    if ax <= 0 or ay <= 0:
        raise ValueError("Lengths ax, ay must be positive.")
    # Separable polynomial (polynomial in x and y for any integer k >= 0)
    px = 1 - ((x - x0)/ax)**k if k > 0 else 1
    py = 1 - ((y - y0)/ay)**k if k > 0 else 1
    p = sp.expand(px * py)

    # Area integral of p over the rectangle
    I = sp.integrate(
            sp.integrate(p, (x, x0 - ax, x0 + ax)),
        (y, y0 - ay, y0 + ay)
    )
    # Final velocity profile u = s * p, with s = Q / I
    s = sp.simplify(Q / I)

    u = sp.simplify(s * p)

elif args.shape == "circular":
    # length = [D] is DIAMETER; domain is disk of radius R centered at (x0,y0)
    if len(args.length) != 1:
        raise ValueError("For circular, --length must be one number: [D] (diameter).")
    D = float(args.length[0])
    if D <= 0:
        raise ValueError("Diameter must be positive.")
    R = D / 2.0

    # To keep a polynomial in x,y, use (x^2+y^2)^m with integer m.
    # If k is even, set m=k/2; if k is odd, use m=ceil(k/2). This yields order 2m.
    # (If you need exactly order k as a polynomial, ensure k is even.)
    m = ceil(max(k, 1) / 2)  # for k=0 -> m=1 gives a flat profile p=1 - (r^2/R^2)^1

    r2 = (x - x0)**2 + (y - y0)**2
    p = sp.expand(1 - (r2**m) / (R**(2*m)))

    # Area integral of p over the disk (use polar analytically: I = π R^2 * m/(m+1))
    I = sp.pi * R**2 * sp.Rational(m, m + 1)

    s = sp.simplify(Q / I)
    u = sp.simplify(s * p)

else:
    raise ValueError("Unknown shape.")

# Output
print("\n=== Inlet velocity profile u(x,y) (SymPy expression) ===\n")
print(u)

# Optional: show a 1D cut along the centerline (rectangular only)
if args.show and args.shape == "rectangular":
    # Plot u(x, y0) across the width
    ux_center = sp.simplify(sp.expand(u.subs(y, y0)))
    print("\nPlotting centerline u(x, y0) with SymPy...")
    sp.plot(ux_center, (x, x0 - ax, x0 + ax),
            title="Centerline x velocity profile",
            xlabel="x", ylabel="u(x, y0)")

    uy_center = sp.simplify(sp.expand(u.subs(x, x0)))
    print("\nPlotting centerline u(x0, y) with SymPy...")
    sp.plot(uy_center, (y, y0 - ay, y0 + ay),
            title="Centerline y velocity profile",
            xlabel="y", ylabel="u(x0, y)")