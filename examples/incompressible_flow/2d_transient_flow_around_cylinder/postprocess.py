import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
import sys

currentPath = sys.argv[1]

filename = currentPath + '/force.00.dat'

f = open(filename, 'r')

lines = f.readlines()[2001:]

data = np.empty((len(lines),4))
for i in range(len(lines)):
    data[i] = lines[i].split()

CD = 2.0*data[:,1]
CL = 2.0*data[:,2]

averagedCD = np.mean(CD)
averagedCL = np.mean(CL)

varCD = np.mean([abs(averagedCD - np.min(CD)), np.max(CD) - averagedCD])
varCL = np.mean([abs(averagedCL - np.min(CL)), np.max(CL) - averagedCL])

print("CD = " + str(averagedCD) + " +/- " + str(varCD))
print("CL = " + str(averagedCL) + " +/- " + str(varCL))

X = 2.0*fft(CL)/len(CL)

freq = fftfreq(len(X), d=0.05)

N = len(freq)

# Get the one-sided specturm
n_oneside = N//2
# get the one side frequency
f_oneside = freq[:n_oneside]

plt.figure(figsize = (12, 6))

plt.loglog(f_oneside, np.abs(X[:n_oneside]), 'b')
plt.xlabel('Freq [Hz]')
plt.ylabel('FFT Amplitude [-]')
plt.grid(which='major')
plt.grid(which='minor', ls=':')
plt.axis([0.01,10,0.00001,1])
plt.savefig(sys.argv[1] + '/cylinderFFT.png')
plt.show()
