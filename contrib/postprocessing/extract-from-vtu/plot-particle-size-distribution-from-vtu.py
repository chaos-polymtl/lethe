#############################################################################
"""
Extraction particles dismeter from a .pvd file and plot the cumulative
distribution function (CDF) based on number using normal distribution.
"""
#############################################################################
'''Importing Libraries'''
import sys
import numpy as np
from lethe_pyvista_tools import *
from scipy.stats import norm
from matplotlib import pyplot as plt
#############################################################################

def cdf_normal(x, mu, sigma):
    """CDF using normal mean and std (mu, sigma"""
    return norm.cdf((x - mu) / sigma)

def cdf_lognormal_X(x, mu_V, sigma_V):
    """CDF using lognormal mean and std (mu_V, sigma_V)"""
    sigma = np.sqrt(np.log(1 + (sigma_V/mu_V)**2))
    mu = np.log(mu_V) - 0.5 * sigma**2
    return norm.cdf((np.log(x) - mu) / sigma)

#############################################################################
# Create the particle object
prm_file_name = "insert_volume_normal_distribution.prm"
pvd_name = 'out.pvd'
particle = lethe_pyvista_tools("./", prm_file_name, pvd_name)
#############################################################################

df = particle.get_df(-1)
mu_i = particle.prm_dict["average diameter"]
sigma_i = particle.prm_dict["standard deviation"]

# Diameter
diameters = df["diameter"][:]

dia_min = diameters.min()
dia_max = diameters.max()

d_plot = np.linspace(dia_min* 0.9, dia_max * 1.1, 200)
N_plot = np.zeros_like(d_plot)

for i,d in enumerate(d_plot):
    N_plot[i] = np.sum(diameters <= d) / len(diameters)

plt.plot(d_plot, N_plot)
plt.plot(d_plot, cdf_normal(d_plot, mu_i, sigma_i), "--")
plt.show()
