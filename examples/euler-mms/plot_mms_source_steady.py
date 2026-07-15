import sympy as sp


def to_lethe(expr):
    """
    Convert SymPy expression to deal.II / Lethe FunctionParser style.
    """
    expr = sp.simplify(expr)
    text = sp.sstr(expr)
    text = text.replace("**", "^")
    return text


def get_steady_mms_expressions():
    x, y, t = sp.symbols("x y t")
    rho_s, beta = sp.symbols("rho_s beta")

    pi = sp.pi

    # ------------------------------------------------------------------
    # Steady manufactured solution
    # ------------------------------------------------------------------
    s = sp.sin(pi * x) * sp.sin(pi * y)

    alpha = sp.Rational(3, 10) + sp.Rational(1, 20) * s
    u = sp.Rational(1, 5) * s
    v = sp.Rational(3, 20) * s

    # Manufactured fluid velocity used in drag
    uf_x = sp.Integer(1)
    uf_y = sp.Integer(0)

    # ------------------------------------------------------------------
    # Solid continuity source
    #
    # d(alpha)/dt + d(alpha*u)/dx + d(alpha*v)/dy = S_alpha
    #
    # Since alpha is steady:
    # d(alpha)/dt = 0
    # ------------------------------------------------------------------
    S_alpha = (
        sp.diff(alpha, t)
        + sp.diff(alpha * u, x)
        + sp.diff(alpha * v, y)
    )

    # ------------------------------------------------------------------
    # Solid x-momentum source
    #
    # rho_s d(alpha*u)/dt
    # + rho_s d(alpha*u*u)/dx
    # + rho_s d(alpha*u*v)/dy
    # + beta alpha (u - uf_x)
    # = S_mx
    #
    # Since alpha and u are steady:
    # d(alpha*u)/dt = 0
    # ------------------------------------------------------------------
    S_mx = (
        rho_s * sp.diff(alpha * u, t)
        + rho_s * sp.diff(alpha * u * u, x)
        + rho_s * sp.diff(alpha * u * v, y)
        + beta * alpha * (u - uf_x)
    )

    # ------------------------------------------------------------------
    # Solid y-momentum source
    #
    # rho_s d(alpha*v)/dt
    # + rho_s d(alpha*v*u)/dx
    # + rho_s d(alpha*v*v)/dy
    # + beta alpha (v - uf_y)
    # = S_my
    #
    # Since alpha and v are steady:
    # d(alpha*v)/dt = 0
    # ------------------------------------------------------------------
    S_my = (
        rho_s * sp.diff(alpha * v, t)
        + rho_s * sp.diff(alpha * v * u, x)
        + rho_s * sp.diff(alpha * v * v, y)
        + beta * alpha * (v - uf_y)
    )

    # Order used by the C++ source-term parser:
    # S_mx ; S_my ; S_alpha
    source_expression = (
        f"{to_lethe(S_mx)}; \n"
        f"{to_lethe(S_my)}; \n"
        f"{to_lethe(S_alpha)}"
    )

    # Order used by analytical solution parser:
    # u ; v ; alpha
    analytical_expression = (
        f"{to_lethe(u)}; \n"
        f"{to_lethe(v)}; \n"
        f"{to_lethe(alpha)}"
    )

    # Initial condition
    initial_u = to_lethe(u)
    initial_v = to_lethe(v)
    initial_alpha = to_lethe(alpha)

    # Boundary values on [-1,1] x [-1,1]:
    # sin(pi*x) or sin(pi*y) is zero on all boundaries
    boundary_u = "0"
    boundary_v = "0"
    boundary_alpha = "0.30"

    return {
        "SOURCE_EXPRESSION": source_expression,
        "ANALYTICAL_EXPRESSION": analytical_expression,
        "INITIAL_U": initial_u,
        "INITIAL_V": initial_v,
        "INITIAL_ALPHA": initial_alpha,
        "BOUNDARY_U": boundary_u,
        "BOUNDARY_V": boundary_v,
        "BOUNDARY_ALPHA": boundary_alpha,
        "RAW_S_ALPHA": sp.simplify(S_alpha),
        "RAW_S_MX": sp.simplify(S_mx),
        "RAW_S_MY": sp.simplify(S_my),
    }


if __name__ == "__main__":
    expressions = get_steady_mms_expressions()

    print("\n--- Source term expression: S_mx ; S_my ; S_alpha ---")
    print(expressions["SOURCE_EXPRESSION"])

    print("\n--- Analytical solution expression: u ; v ; alpha ---")
    print(expressions["ANALYTICAL_EXPRESSION"])

    print("\n--- Initial conditions ---")
    print("u =", expressions["INITIAL_U"])
    print("v =", expressions["INITIAL_V"])
    print("alpha =", expressions["INITIAL_ALPHA"])

    print("\n--- Boundary conditions ---")
    print("u =", expressions["BOUNDARY_U"])
    print("v =", expressions["BOUNDARY_V"])
    print("alpha =", expressions["BOUNDARY_ALPHA"])