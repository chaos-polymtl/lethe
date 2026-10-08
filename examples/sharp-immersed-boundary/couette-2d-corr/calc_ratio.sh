
import math


def get_fx(filename):
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    header = lines[0].split()
    col = header.index("f_x")
    last_line = lines[-1].split()

    return float(last_line[col])


Fx0 = get_fx("out_ref/force_ref.02.dat")
Fx = get_fx("out/force_7p_eps0.02.dat")

print("========== FORCES ==========")
print(f"Fx0 (sans particules) = {Fx0:.6f}")
print(f"Fx  (avec particules) = {Fx:.6f}")


mu_fluide = 10.0
rho_fluide = 1.0
rho_particule = 1.0

N_particules = 7
r_particule = 1.0

Lx = 50.0
Ly = 20.0


A_domaine = Lx * Ly
A_particules = N_particules * math.pi * r_particule**2

phi = A_particules / A_domaine

print("\n========== FRACTION DE PARTICULES ==========")
print(f"A domaine       = {A_domaine:.6f}")
print(f"A particules    = {A_particules:.6f}")
print(f"phi             = {phi:.6f}")
print(f"phi (%)         = {phi * 100:.3f} %")


mu_einstein = mu_fluide * (1.0 + 2.5 * phi)

augmentation_einstein = (
    (mu_einstein - mu_fluide) / mu_fluide
) * 100.0

print("\n========== EINSTEIN ==========")
print(f"mu fluide       = {mu_fluide:.6f}")
print(f"mu Einstein     = {mu_einstein:.6f}")
print(f"augmentation    = {augmentation_einstein:.3f} %")


ratio_force = Fx / Fx0

print("\n========== SIMULATION ==========")
print(f"Fx / Fx0        = {ratio_force:.6f}")

mu_sim = mu_fluide * ratio_force

augmentation_sim = (
    (mu_sim - mu_fluide) / mu_fluide
) * 100.0

print(f"mu simulation   = {mu_sim:.6f}")
print(f"augmentation    = {augmentation_sim:.3f} %")


difference = mu_sim - mu_einstein

erreur_relative = (
    abs(mu_sim - mu_einstein) / mu_einstein
) * 100.0

print("\n========== COMPARAISON ==========")
print(f"mu Einstein     = {mu_einstein:.6f}")
print(f"mu simulation   = {mu_sim:.6f}")
print(f"difference      = {difference:.6f}")
print(f"erreur relative = {erreur_relative:.3f} %")


rho_eff = (
    (1.0 - phi) * rho_fluide
    + phi * rho_particule
)

print("\n========== MASSE VOLUMIQUE ==========")
print(f"rho effective   = {rho_eff:.6f}")
