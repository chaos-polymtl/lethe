import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import argparse

colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.markersize'] = '11'
plt.rcParams['markers.fillstyle'] = "none"
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.columnspacing'] = 2
plt.rcParams['legend.handlelength'] = 3
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = False
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['font.size'] = '25'
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['font.serif']='cm'
plt.rcParams['savefig.bbox']='tight'

parser = argparse.ArgumentParser(description='Arguments for the post-processing')
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)

B = 2
vol = B**3
rho = 0.00759
rho_0= 0.00124
mass = vol * rho
I= mass*B**2 / 6
g=981

Uc = np.sqrt(4./3. * g * B * np.abs(rho-rho_0)/rho_0 )
tc = B / Uc 
t_fact = 1/tc


mat=np.loadtxt("out/ib_force.00.dat",skiprows=1)
t=mat[:,1]
vy=np.abs(mat[:,15])

experimental_data = np.loadtxt("S18.dat",skiprows=1)
plt.plot(experimental_data[:,0],experimental_data[:,1],'ks', mfc='none',label="Wang et al.")
plt.plot(t*t_fact,vy,color=colors[0],lw=3, label="Lethe")

plt.ylabel("Sedimentation velocity, $v_y$ [cm/s]")
plt.xlabel("Dimensionless time, $t^*$ [-]")
plt.legend()
if (parser.validate):
    plt.savefig("cuboid-sedimentation-velocity.pdf")
    plt.close()

    lethe_validation_data = np.column_stack((t*t_fact, vy))
    np.savetxt("cuboid-sedimentation-velocity.dat", lethe_validation_data, header="t* vy")
else:
    plt.savefig("velocity-comparison.png", dpi=300)
    plt.show()