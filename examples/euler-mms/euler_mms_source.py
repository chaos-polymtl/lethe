import sympy as sp


def to_lethe(expr):
    """
    Convert SymPy expression to deal.II / Lethe FunctionParser style.
    """
    expr = sp.simplify(expr)
    text = sp.sstr(expr)
    text = text.replace("**", "^")
    return text


def get_mms_expressions():
    x, y, t = sp.symbols("x y t")
    rho_s, beta = sp.symbols("rho_s beta")

    pi = sp.pi

    # ------------------------------------------------------------------
    # Manufactured solution
    # ------------------------------------------------------------------
    alpha = (
        sp.Rational(3, 10)
        + sp.Rational(1, 20)
        * sp.sin(pi * x)
        * sp.sin(pi * y)
        * sp.sin(2 * pi * t)
    )

    u = (
        sp.Rational(1, 5)
        * sp.sin(pi * x)
        * sp.sin(pi * y)
        * sp.cos(2 * pi * t)
    )

    v = (
        sp.Rational(3, 20)
        * sp.sin(pi * x)
        * sp.sin(pi * y)
        * sp.sin(2 * pi * t)
    )

    # Manufactured fluid velocity used in drag
    uf_x = sp.Integer(1)
    uf_y = sp.Integer(0)

    # ------------------------------------------------------------------
    # Solid continuity source
    #
    # d(alpha)/dt + d(alpha*u)/dx + d(alpha*v)/dy = S_alpha
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
        f"{to_lethe(S_mx)}; \n "
        f"{to_lethe(S_my)}; \n "
        f"{to_lethe(S_alpha)}"
    )

    # Order used by analytical solution parser:
    # u ; v ; alpha
    analytical_expression = (
        f"{to_lethe(u)}; \n " 
        f"{to_lethe(v)}; \n " 
        f"{to_lethe(alpha)}"
    )

    # Initial condition at t = 0
    initial_u = to_lethe(u.subs(t, 0))
    initial_v = to_lethe(v.subs(t, 0))
    initial_alpha = to_lethe(alpha.subs(t, 0))

    # Boundary values on unit square:
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
    }


if __name__ == "__main__":
    expressions = get_mms_expressions()

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