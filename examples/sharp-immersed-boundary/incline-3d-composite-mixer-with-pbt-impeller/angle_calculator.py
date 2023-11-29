#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 15:35:32 2021
"""

#!/usr/bin/python

from sympy import *
import numpy as np


x, y, z,  ox, oy, oz, pi, t= symbols('x y z ox oy oz pi t ')

def rot_axisx(theta):
    """Returns a rotation matrix for a rotation of theta (in radians) about
    the 1-axis.
    [...]
    """
    ct = cos(theta)
    st = sin(theta)
    lil = ((1, 0, 0),
           (0, ct, -st),
           (0, st, ct))
    return Matrix(lil)


def rot_axisy(theta):
    """Returns a rotation matrix for a rotation of theta (in radians) about
    the 2-axis.
    [...]
    """
    ct = cos(theta)
    st = sin(theta)
    lil = ((ct,0,st),
           (0, 1, 0),
           (-st, 0, ct))
    return Matrix(lil)

def rot_axisz(theta):
    """Returns a rotation matrix for a rotation of theta (in radians) about
    the 3-axis.
    [...]
    """
    ct = cos(theta)
    st = sin(theta)
    lil = ((ct, -st, 0),
           (st, ct, 0),
           (0, 0, 1))
    return Matrix(lil)


def rotation_matrix_to_xyz_angles(R):
    """
    Extracts XYZ rotation angles from a given rotation matrix.

    Parameters:
    R (Matrix): A 3x3 rotation matrix.

    Returns:
    tuple: A tuple of rotation angles (theta_x, theta_y, theta_z) in radians.
    """
    if R.shape != (3, 3):
        raise ValueError("Input must be a 3x3 matrix.")

    # Calculating the angles
    theta_x = atan2(-R[1, 2], R[2, 2])
    theta_y = asin(R[0, 2])
    theta_z = atan2(-R[0, 1], R[0, 0])

    return theta_x, theta_y, theta_z


# Rotation matrix for a small time step dt
initial_rot_x=0
initial_rot_y=pi/4
initial_rot_z=0

Initial_rotation=rot_axisx(initial_rot_x)*rot_axisy(initial_rot_y)*rot_axisz(initial_rot_z)

# Angular velocity vector
ox=-1*np.pi*2*np.sqrt(2)/2.0
oy=0
oz=-1*np.pi*2*np.sqrt(2)/2.0

# Magnitude of the angular velocity vector
omega_mag = sqrt(ox**2 + oy**2 + oz**2)

# Unit vector along the direction of angular velocity
u_x = ox / omega_mag
u_y = oy / omega_mag
u_z = oz / omega_mag

# Rodrigues' rotation formula components
K = Matrix([[0, -u_z, u_y],
            [u_z, 0, -u_x],
            [-u_y, u_x, 0]])

I = Matrix([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])


R = I + sin(omega_mag*t) * K + (1 - cos(omega_mag*t)) * K**2


theta_x, theta_y, theta_z=rotation_matrix_to_xyz_angles(R*Initial_rotation)

# Print orientation
print(str(theta_x).replace("**","^")+';'+str(theta_y).replace("**","^")+';'+str(theta_z).replace("**","^"))