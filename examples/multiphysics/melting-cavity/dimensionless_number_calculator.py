import numpy as np

#Target dimensionless numbers
Ra=1e5 #Rayleigh Number
Gr=5.5e6
St=0.041

H=0.714
rho=1
Tw=38
Ts=29.8
g=1
beta=1
cp=1

# We fix enthalpy using Stefan number
hl = cp * (Tw-Ts)/St
print("hl value is: ", hl)

# We fix viscosity to get Grashoff
mu = np.sqrt(g*beta*(Tw-Ts)*H**3/Gr)
print("dynamic viscosity is: ", mu)

# We fix conductivity to get Rayleigh
k = H**3 * (Tw-Ts) * g * beta / mu / Ra
print("Conductivity should be: ", k)


# The final dimensionless numbers are:
calc_Ra= H**3 * (Tw-Ts) * g * beta / mu / k
calc_Gr= g*beta*(Tw-Ts)*H**3/mu**2
calc_St=cp * (Tw-Ts) / hl

print("Calc Ra: ", calc_Ra)
print("Calc Gr: ", calc_Gr)
print("Calc St: ", calc_St)
print("Final time: ", 0.1/k/St)

