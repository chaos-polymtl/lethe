import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import optimize
import jinja2
import os

## This small code provides the Young's modulus of the particle to optain a certain value of normal rigidity constant.
## To run this code correcly, you just need to enter the next three sets of inputs and run it.

class Wall_prop:
    poisson_coef_wall = 0.3    # Wall poisson ratio
    young_wall = 1000000000000 # Wall Young's modulus


class Particle_prop:
    poisson_coef_particle = 0.3                       # Particle poisson ratio
    radius_particle = 0.1                             # Particle radius
    density = 2600                                    # Particle density
    volume = (4 / 3) * np.pi * radius_particle  ** 3  # Particle volume
    mass_particle = volume * density                  # Particle mass


class Param:
    tol = 1
    N_max = 10000
    kn_wanted = float(sys.argv[1])  # 5 * (10**5)   //   5 * (10**4)
    y_estim_1 = 1e3          # a , Lowerbound
    y_estim_2 = 1e10         # b , Upperbound


def Young_effective(Poisson_wall, Young_wall, Poisson_particle, Young_particle):
    # Calculate the effective Young's modulus between two materials
    Y_effective = ((1-Poisson_wall**2)/Young_wall + (1-Poisson_particle**2)/Young_particle)**(-1)
    return Y_effective


def f(x):
    y_eff = Young_effective(Wall_prop.poisson_coef_wall,Wall_prop.young_wall,Particle_prop.poisson_coef_particle,x)
    parentheses = (15 * Particle_prop.mass_particle * 1.**2. / (16 * Particle_prop.radius_particle**0.5 * y_eff ))**0.2
    return (16/15) * Particle_prop.radius_particle**0.5 * y_eff * parentheses - Param.kn_wanted


def copy_and_replace_template(source_file, destination_file, ER, er, K):
    with open(source_file, 'r') as f:
        content = f.read()

print("The young's modulus of the particle must be equal to : ")
yp = optimize.bisect(f,Param.y_estim_1,Param.y_estim_2)           # Young's modulus of the particle
print(yp)

# Figure
Young_particle = np.linspace(1e4,1e10,1000)
kn1 = []
kn2 = []
for i in Young_particle:
    kn2 = np.append(kn2,f(i)+Param.kn_wanted)

plt.loglog(Young_particle, kn2)
plt.loglog(yp, Param.kn_wanted,"sr")
plt.grid()
plt.title("Normal stiffness constant as a function\n of the Young's modulus of the particle")
plt.xlabel("Young's modulus of the particle")
plt.ylabel("Normal stiffness constant")
plt.show()


# Case generator
PATH = os.getcwd()
CASE_PREFIX = 'bouncing_particle_restitution-coef_'
PRM_FILE = 'bouncing_particle_original.prm'
TEMPLATE_FOLDER = 'template'

# System characteristics
er = np.array([0.5,0.6,0.7,0.8,0.9,1.0])      # Restitution coefficient

templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(PRM_FILE)

for val in er:

    output_text = template.render(ER=f"{int(val*10):02d}", er=val, YP=yp)
    prm_file_name = f"bouncing_particle_{int(val*10):02d}.prm"

    # Write the output text to the prm file
    output_file_path = os.path.join("./", prm_file_name)
    with open(output_file_path, 'w') as f:
        f.write(output_text)

