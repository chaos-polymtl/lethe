import numpy as np
import matplotlib.pyplot as plt

D = 0.27
N = 10 / 2 / np.pi

Re = np.logspace(-1,np.log10(100),25)

# Simulation Data
C,Tx,Ty,Tz = np.loadtxt("gather.dat",unpack=True)
n_dat = np.size(Tz)
D = 0.27
Np = 2 * np.pi * np.abs(Tz) / 1 / N**2 / D**5

# Experimental data
Re_1,Np_1, Re_2,Np_2, Re_3,Np_3 = np.loadtxt("experimental.dat",unpack=True,skiprows=1)

plt.loglog(Re_1,Np_1,'o',mfc='None',label="Experiment 1")
plt.loglog(Re_2,Np_2,'*',mfc='None',label="Experiment 2")
plt.loglog(Re_3,Np_3,'^',mfc='None',label="Experiment 3")
plt.loglog(Re[:n_dat],Np,'-',label="Lethe",color="black")

plt.xlabel("$Re$")
plt.ylabel("$N_p$")
plt.legend()

plt.show()
