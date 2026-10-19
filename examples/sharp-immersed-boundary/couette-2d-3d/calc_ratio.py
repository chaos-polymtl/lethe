#Compare wall force (and effective viscosity) with and without particles
#against the Einstein dilute-suspension model in 2D Couette flow.

import argparse
import math


def get_fx(filename):
    """Read the last data row of a Lethe force file and return f_x."""
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    header = lines[0].split()
    col = header.index("f_x")
    last_line = lines[-1].split()
    return float(last_line[col])


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare Couette wall forces with/without particles "
            "to the Einstein viscosity model."
        )
    )
    parser.add_argument(
        "force_ref",
        help="Force file without particles (e.g. out_ref/force_ref.02.dat)",
    )
    parser.add_argument(
        "force_particles",
        help="Force file with particles (e.g. out_2d/force_7p_eps0.02.dat)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # --- Wall forces from Lethe output ---------------------------------
    Fx0 = get_fx(args.force_ref)          # empty channel (reference)
    Fx = get_fx(args.force_particles)     # channel with fixed particles

    print("========== FORCES ==========")
    print(f"Fx0 (no particles)   = {Fx0:.6f}")
    print(f"Fx  (with particles) = {Fx:.6f}")

    # --- Physical and geometric parameters -----------------------------
    mu_fluid = 10.0
    rho_fluid = 1.0
    rho_particle = 1.0

    n_particles = 7
    r_particle = 1.0

    Lx = 50.0
    Ly = 20.0

    # Area fraction of particles in the 2D domain
    A_domain = Lx * Ly
    A_particles = n_particles * math.pi * r_particle**2
    phi = A_particles / A_domain

    print("\n========== PARTICLE FRACTION ==========")
    print(f"Domain area     = {A_domain:.6f}")
    print(f"Particle area   = {A_particles:.6f}")
    print(f"phi             = {phi:.6f}")
    print(f"phi (%)         = {phi * 100:.3f} %")

    # --- Einstein dilute-limit viscosity -------------------------------
    # mu_E = mu (1 + 5/2 phi)  for spheres in 3D; used here as reference
    mu_einstein = mu_fluid * (1.0 + 2.5 * phi)
    increase_einstein = ((mu_einstein - mu_fluid) / mu_fluid) * 100.0

    print("\n========== EINSTEIN ==========")
    print(f"mu fluid        = {mu_fluid:.6f}")
    print(f"mu Einstein     = {mu_einstein:.6f}")
    print(f"increase        = {increase_einstein:.3f} %")

    # --- Effective viscosity inferred from wall force ------------------
    # In Couette flow, wall shear force scales with viscosity, so
    # mu_sim / mu_fluid = Fx / Fx0
    force_ratio = Fx / Fx0
    mu_sim = mu_fluid * force_ratio
    increase_sim = ((mu_sim - mu_fluid) / mu_fluid) * 100.0

    print("\n========== SIMULATION ==========")
    print(f"Fx / Fx0        = {force_ratio:.6f}")
    print(f"mu simulation   = {mu_sim:.6f}")
    print(f"increase        = {increase_sim:.3f} %")

    # --- Comparison ----------------------------------------------------
    difference = mu_sim - mu_einstein
    relative_error = abs(mu_sim - mu_einstein) / mu_einstein * 100.0

    print("\n========== COMPARISON ==========")
    print(f"mu Einstein     = {mu_einstein:.6f}")
    print(f"mu simulation   = {mu_sim:.6f}")
    print(f"difference      = {difference:.6f}")
    print(f"relative error  = {relative_error:.3f} %")

    # Mixture density (particles and fluid have equal density here)
    rho_eff = (1.0 - phi) * rho_fluid + phi * rho_particle

    print("\n========== EFFECTIVE DENSITY ==========")
    print(f"rho effective   = {rho_eff:.6f}")


if __name__ == "__main__":
    main()
