# %%
# pip install sympy

# %%
from sympy import *

x, y, z, nu = symbols('x y z nu')

# Viscous term
def deriv2(f):
    d2dx = diff(f,x,x)
    d2dy = diff(f,y,y)

    return d2dx+d2dy

# Advection term
def convec(f,u,v):
    cx = u*diff(f,x)
    cy = v*diff(f,y)

    return cx+cy
 
# Pressure gradient (the pressure is P/\rho)
def dx(f):
    return diff(f,x)

def dy(f):
    return diff(f,y)

u = -2*(sin(pi*x))**2*sin(pi*y)*cos(pi*y)
v = 2*sin(pi*x)*cos(pi*x)*(sin(pi*y))**2

p = sin(pi*x)*sin(pi*y)

# Calculating the source term components

# x-component of the momentum equation: this is the first term in the source terms of the fluid dynamics subsection. See https://chaos-polymtl.github.io/lethe/documentation/parameters/cfd/source_term.html#source-term
Q_x = convec(u,u,v) + dx(p) - nu*(deriv2(u))
print("x-component of momentum conservation equation source term:", simplify(Q_x))
# y-component of the momentum equation: this is the second term in the source terms of the fluid dynamics subsection
Q_y = convec(v,u,v) + dy(p) - nu*(deriv2(v))
print("x-component of momentum conservation equation source term:", simplify(Q_y))
# Continuity equation
Q_c = dx(u)+dy(v)
print("Continuity equation source term:", Q_c)

# Since we are calculating the source term that corresponds to whichever solution we are imposing, the typical Poisson equation is no longer valid, so this is useless
# print("Laplacian of the pressure")
# print(simplify(deriv2(p)))
# print("Velocity term in the Poisson equation")
# velo_poisson=simplify(dx(convec(u,u,v))+dy(convec(v,u,v)))
# print(velo_poisson)
# velo_poisson=trigsimp(dx(convec(u,u,v))+dy(convec(v,u,v)))
# print("With trigsimp")
# print(velo_poisson)
# print(u,v, p)


