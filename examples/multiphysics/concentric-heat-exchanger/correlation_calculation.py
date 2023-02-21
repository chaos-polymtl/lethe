import numpy as np
import matplotlib.pyplot as plt

def area(r):
    return r*r*np.pi

u_in = 1e-3 # m/s
u_out = 0.2e-3 # m/s

R_in = 0.002
R_ext = 0.003
R_out = 0.005
L = 1

# Calculate inner and outer flow rates
A_in = area(R_in)
A_ext = area(R_ext)
A_out = area(R_out)

Q_in = u_in * A_in
Q_out = u_out * A_out

print("Respective flow rates are (Q_in): ", Q_in, " (Q_out): ", Q_out )

# Calculate h coefficient. Assume uniform q'' (Nu_d=4.36)

k_fluid=5
h_in  = k_fluid *4.36 / 2 / R_in
h_out = k_fluid *4.36 / 2 / R_out

# Calculate conductivity coefficient
k_solid=45
h_tube =  k_solid  / np.log(R_ext/R_in) / R_in

print("Respective h coefficients  (h_in): ", h_in, " (h_out): ", h_out, " (h_tube): ", h_tube )

# Calculate overall coefficient at inner surface
U = 1 / (1/h_in + 1/h_out + 1/h_tube)

print ("Global heat transfer coefficient U: ", U, " W/m^2")
print ("Global heat transfer coefficient U: ", U* R_in*2*L, " W ")



