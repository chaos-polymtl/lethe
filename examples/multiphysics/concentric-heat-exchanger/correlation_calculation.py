import numpy as np
import matplotlib.pyplot as plt

def area(r1,r2):
    return r2*r2*np.pi-r1*r1*np.pi

u_in = 10e-3 # m/s
u_out = 4e-3 # m/s
Cp = 4180
rho = 1000
k_fluid=0.6
k_solid=398
mu_fluid=1e-6
Pr=6.9

T_hi=100
T_ci=0

R_in = 0.001
R_ext = 0.002
R_out = 0.003
L = 0.05

# Calculate inner and outer flow rates
A_in = area(0,R_in)
A_ext = area(R_in,R_ext)
A_out = area(R_ext,R_out)

Q_in = u_in * A_in
Q_out = u_out * A_out

print("Respective flow rates are (Q_in): ", Q_in, " (Q_out): ", Q_out )

# Calculate h coefficient. Assume uniform q'' (Nu_d=4.36)

Re_in = rho * u_in * 2 * R_in / mu_fluid
Re_out = rho * u_out * 2 * (R_out-R_ext) / mu_fluid

Nu_in = 1.86 * (Re_in * Pr /L * 2 * R_in)**0.33
Nu_out = 1.86 * (Re_in * Pr /L * 2 * (R_out-R_ext))**0.33

h_in  = k_fluid *Nu_in / 2 / R_in
h_out = k_fluid *Nu_out / 2 / R_out

# Calculate conductivity coefficient
h_tube =  k_solid  / np.log(R_ext/R_in) / R_in

print("Respective h coefficients  (h_in): ", h_in, " (h_out): ", h_out, " (h_tube): ", h_tube )

# Calculate overall coefficient at inner surface
U = 1 / (1/h_in + 1/h_out + 1/h_tube)

print ("Global heat transfer coefficient U: ", U, " W/m^2/Celsius")
print ("Global heat transfer coefficient U: ", U* R_in*2*L, " W/Celsius ")

#Use NTU method to find temperature
C_in = Q_in * Cp * rho
C_out = Q_out * Cp * rho
Cmin = np.min([C_in,C_out])
Cmax = np.max([C_in,C_out])
Cr = Cmin/Cmax
print("Cr is: ", Cr)
NTU = U * R_in*2*L / Cmin

print("NTU is: ", NTU)
eps = (1 - np.exp(-NTU*(1-Cr))) / (1-Cr*np.exp(-NTU*(1-Cr)))

print("epsilon is: ", eps)
T_ho = T_hi - eps * (T_hi-T_ci)
print("Outlet temperature should be: ", T_ho )

print("Enthalpies on hot side should be in: ", C_in * T_hi, " out: ", C_in * T_ho)







