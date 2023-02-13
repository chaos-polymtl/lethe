#############################################################################
"""
Postprocessing code for rising-bubble example

"""
#############################################################################

'''Importing Libraries'''
import numpy as np
import sys
import matplotlib.pyplot as plt
#############################################################################

#Take case path as argument and store it
filename = sys.argv[1]
t,x,y,vx,vy=np.loadtxt(filename,skiprows=1,unpack=True)




#Data from Zahedi, Kronbichler and Kreiss (2012)
x_ref = [0.047, 0.139, 0.232 ,0.324 ,0.416 ,0.508 ,0.601 ,0.693 ,0.785 ,0.878 ,0.97 ,1.062 ,1.154 ,1.247 ,1.339 ,1.431 ,1.524 ,1.616 ,1.708 ,1.8 ,1.893 ,1.985 ,2.077 ,2.17 ,2.262 ,2.354 ,2.446 ,2.539 ,2.631 ,2.723 ,2.816 ,2.908 ,2.979]
y_ref = [0.501 ,0.505 ,0.513 ,0.525 ,0.54 ,0.557 ,0.576 ,0.597 ,0.618 ,0.641 ,0.663 ,0.685 ,0.707 ,0.729 ,0.75 ,0.77 ,0.791 ,0.811 ,0.83 ,0.849 ,0.868 ,0.886 ,0.904 ,0.922 ,0.94 ,0.958 ,0.976 ,0.993 ,1.011 ,1.029 ,1.047 ,1.064 ,1.078]
x_vel = [0.001 ,0.03 ,0.056 ,0.081 ,0.106 ,0.136 ,0.17 ,0.203 ,0.237 ,0.271 ,0.304 ,0.338 ,0.372 ,0.41 ,0.452 ,0.494 ,0.544 ,0.607 ,0.687 ,0.779 ,0.872 ,0.964 ,1.056 ,1.149 ,1.241 ,1.333 ,1.426 ,1.518 ,1.61 ,1.702 ,1.795 ,1.887 ,1.979 ,2.072 ,2.164 ,2.256 ,2.348 ,2.441 ,2.533 ,2.625 ,2.718 ,2.81 ,2.902 ,2.978]
y_vel = [0.004 ,0.017 ,0.029 ,0.041 ,0.053 ,0.066 ,0.081 ,0.095 ,0.109 ,0.122 ,0.135 ,0.147 ,0.158 ,0.17 ,0.182 ,0.194 ,0.205 ,0.218 ,0.229 ,0.237 ,0.24 ,0.241 ,0.239 ,0.236 ,0.231 ,0.227 ,0.222 ,0.217 ,0.213 ,0.208 ,0.205 ,0.201 ,0.198 ,0.196 ,0.194 ,0.192 ,0.192 ,0.191 ,0.191 ,0.192 ,0.192 ,0.193 ,0.193 ,0.194]

fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
ax0.plot(t, y, '-k', label="Simulation")
ax0.plot(x_ref, y_ref, 'ro',label="Reference - Zahedi, Kronbichler and Kreiss (2012)")
ax0.set_ylabel(r'Bubble center height')
ax0.set_xlabel(r'$t$')
ax0.legend(loc="upper left")
fig0.savefig(f'./ymean-t.png')
plt.show()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(t, vy , '-k', label="Simulation")
ax1.plot(x_vel, y_vel, 'ro',label="Reference - Zahedi, Kronbichler and Kreiss (2012)")
ax1.set_ylabel(r'Rise velocity')
ax1.set_xlabel(r'$t$')
ax1.legend(loc="upper left")
ax1.legend(loc=4)
fig1.savefig(f'./bubble-rise-velocity.png')
plt.show()
