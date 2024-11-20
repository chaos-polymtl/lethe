# SPDX-FileCopyrightText: Copyright (c) 2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq,fftshift
from numpy.fft import fft,ifft


#Used in read_my_data function
def is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    return True

#Function to smooth a noisy signal based on the Savitzky-Golay filter. All details are explained in the function

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(int(window_size))
        order = np.abs(int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")

    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.asmatrix([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')

#Function allowing to turn .dat files outputed by the simulations into readable python arrays, Thanks to Lucka Barbeau
# for saving me the trouble of doing this work !
def read_my_data(results_path):
    force_file = open(results_path, 'r')
    list_of_list_of_vars_name = [[]];
    list_of_list_of_vars = [[]];
    fix_vars_name = False

    nb_set_of_vars = 0;
    for line in force_file.readlines():
        i = 0;
        for word in line.split():
            if is_number(word):
                fix_vars_name = True
                list_of_list_of_vars[nb_set_of_vars][i] = np.append(list_of_list_of_vars[nb_set_of_vars][i],
                                                                    float(word))
                i += 1
            else:
                if word != 'particle':
                    if fix_vars_name:
                        nb_set_of_vars += 1
                        list_of_list_of_vars_name.append([])
                        list_of_list_of_vars.append([])
                        fix_vars_name = False
                    list_of_list_of_vars_name[nb_set_of_vars].append(word)
                    list_of_list_of_vars[nb_set_of_vars].append(np.array([]))
                    # locals() [word]=np.array([])
    return list_of_list_of_vars_name, list_of_list_of_vars


def printClCd():

    pathList = ["./simu_batch/naca_0.00/force.00.dat",
                "./simu_batch/naca_1.00/force.00.dat",
                "./simu_batch/naca_3.00/force.00.dat",
                "./simu_batch/naca_5.00/force.00.dat",
                "./simu_batch/naca_7.00/force.00.dat",
                "./simu_batch/naca_9.00/force.00.dat",
                "./simu_batch/naca_11.00/force.00.dat",
                "./simu_batch/naca_13.00/force.00.dat",
                "./simu_batch/naca_15.00/force.00.dat"]

    angleList = [0, 1, 3, 5, 7, 9, 11, 13,15]

    rmsCl = []
    rmsCd = []

    meanCl = []
    meanCd = []

    figureCl, axisCl = plt.subplots(5, 2)
    axisCl = axisCl.flatten()
    
    cut = 5000
    cutLast = 200

    for i in range(0,len(pathList)):
        list_of_list_of_vars_name, list_of_list_of_vars = read_my_data(pathList[i])
        list_var_interet = [2 * list_of_list_of_vars[0][1][cut:-cutLast], 2 * list_of_list_of_vars[0][2][cut:-cutLast],
                            list_of_list_of_vars[0][4][cut:-cutLast],
                            list_of_list_of_vars[0][5][cut:-cutLast], list_of_list_of_vars[0][7][cut:-cutLast],
                            list_of_list_of_vars[0][8][cut:-cutLast]]
        for j in range(0, len(list_var_interet)):
            list_var_interet[j] = savitzky_golay(list_var_interet[j], 11, 1)


        rmsCl.append(np.sqrt(np.mean(list_var_interet[1]**2)))
        rmsCd.append(np.sqrt(np.mean(list_var_interet[0]**2)))
        meanCl.append(np.mean(list_var_interet[1]))
        meanCd.append(np.mean(list_var_interet[0]))

        axisCl[i].plot(list_of_list_of_vars[0][0][cut:-cutLast],list_var_interet[1])
        axisCl[i].set_title(str(angleList[i]) + "degrees" + "//Cl_rms = " + f"{np.sqrt(np.mean(list_var_interet[1]**2)):.2e}")

    # Hide the last subplot
    axisCl[-1].axis('off')

    #The following data come from Kouser et al. (2021)

    RefAngles =  [0, 5.052966101694915, 7.023305084745762, 8.040254237288135, 9.025423728813559, 10.042372881355933, 13.029661016949152, 15, 15.985169491525424]
    RefCl = [0, 0.24740532959326786, 0.302945301542777, 0.31977559607293127, 0.3702664796633941, 0.4140252454417952, 0.5570827489481065, 0.6768392370572206, 0.6631136044880784]
    RefCd = [0.1211781206171108, 0.12959326788218792, 0.13800841514726506, 0.14137447405329592, 0.15483870967741933, 0.1649368863955119, 0.2187938288920056, 0.2844686648501362,0.3188010899182561]


    rmsCl = np.array(rmsCl)
    rmsCd = np.array(rmsCd)

    figure,axis = plt.subplots(1,2)
    axis[0].plot(angleList,rmsCl,marker="o",label="Simulated lift coefficient")
    #axis[0].plot(angleList, meanCl, marker="o", label="Numerical C_L") #Uncomment to plot the mean Cl
    axis[0].plot(RefAngles[:-1],RefCl[:-1],marker="*",label="Kouser et. al")
    axis[0].legend(loc='upper left',fontsize=20)
    axis[0].set_title("Lift coefficient, Re = 1000",fontsize=20)
    axis[0].set_xlabel("Angle of attack (°)",fontsize=20)
    axis[0].set_ylabel("RMS value of lift coefficient",fontsize=20)
    axis[0].tick_params(axis='y', labelsize=20)
    axis[0].tick_params(axis='x', labelsize=20)

    axis[1].plot(angleList,rmsCd,marker="o",label="Simulated drag coefficient")
    #axis[1].plot(angleList, meanCd, marker="o", label="Simulated drag coefficient") #Uncomment to plot the mean Cd
    axis[1].plot(RefAngles[:-1], RefCd[:-1], marker="*", label="Kouser et. al")
    #axis[1].plot(angleList, meanCd, marker="o")
    axis[1].legend(loc='upper left',fontsize=20)
    axis[1].set_title("Drag coefficient, Re = 1000",fontsize=20)
    axis[1].set_xlabel("Angle of attack (°)",fontsize=20)
    axis[1].set_ylabel("RMS value of drag coefficient",fontsize=20)
    axis[1].tick_params(axis='y', labelsize=20)
    axis[1].tick_params(axis='x', labelsize=20)
    axis[1].yaxis.set_label_position("right")
    axis[1].yaxis.tick_right()

    plt.show()

printClCd()


#Used for plotting the fft of the lift coefficient for the alpha = 15° case. To perform fft on the other angles of attack,
#simply replace angleIndex =8 by the index of the desired angle.
def spectralAnalysis():

    # Data from Figure 7. of doi/10.1177/17568293211055656
    x_ref = [0.4793106605266697, 0.54137946819437, 0.5879316658809266, 0.6344830743197738, 0.6551724137931041, 0.6758625425141431, 0.6913799417429962, 0.7120692812163251, 0.7224139509529903, 0.7379313501818434, 0.7534487494106952, 0.7844827586206892, 0.8362068965517243, 0.9086211632037985, 0.9706899708714973, 1.0224141088025327, 1.037931508031385, 1.0637935769969018, 1.1000003156990847, 1.1931039218244885, 1.275862068965516, 1.3482763356175909, 1.3689656750909214, 1.3793103448275865, 1.394827744056438, 1.4103451432852914, 1.425862542514143, 1.4775866804451783, 1.5448282176050636, 1.6224144245016157]
    y_ref = [0.002672605415990023, 0.0034744053921120597, 0.004810708100107052, 0.008285093101857979, 0.012828524348077122, 0.025389771842266352, 0.10048998505110408, 0.05077952329417159, 0.020044561010286305, 0.012828524348077122, 0.008285093101857979, 0.0056124876858679555, 0.003207148928585284, 0.0021381026841170293, 0.0021381026841170293, 0.0037416720508193927, 0.010423185590794469, 0.003207148928585284, 0.0008017999761220366, 0, 0.00026727685388790905, 0.0018708462205902724, 0.0056124876858679555, 0.012294001225843013, 0.02458799225650547, 0.006681533930336192, 0.0034744053921120597, 0.0016035795618829206, 0.0010690666348293694, 0.0005345435125952609]

    pathList = ["./simu_batch/naca_0.00/force.00.dat",
                "./simu_batch/naca_1.00/force.00.dat",
                "./simu_batch/naca_3.00/force.00.dat",
                "./simu_batch/naca_5.00/force.00.dat",
                "./simu_batch/naca_7.00/force.00.dat",
                "./simu_batch/naca_9.00/force.00.dat",
                "./simu_batch/naca_11.00/force.00.dat",
                "./simu_batch/naca_13.00/force.00.dat",
                "./simu_batch/naca_15.00/force.00.dat"]

    angleList = [0, 1, 3, 5, 7, 9, 11, 13,15]

    angleIndex = 8

    #The following "cut" parameter allows us to truncate the signal from the beginning to avoid having to deal with unphysical
    #values of drag/lift due to numercal artifacts
    cut = 500;


    #for i in range(0,9):
    list_of_list_of_vars_name, list_of_list_of_vars = read_my_data(pathList[angleIndex])
    list_var_interet = [2 * list_of_list_of_vars[0][1][cut:], 2 * list_of_list_of_vars[0][2][cut:],
                            list_of_list_of_vars[0][4][cut:],
                            list_of_list_of_vars[0][5][cut:], list_of_list_of_vars[0][7][cut:],
                            list_of_list_of_vars[0][8][cut:]]
                            
    #print(list_of_list_of_vars[0][2][500:])  #
    #print(list_of_list_of_vars_name)

    # Sampling rate : since the timestep is non-uniform, choice was made to use a mean timestep and perform
    #a regular fft on the signal. for the sake of simplicity. A decicated signal processing tool exists for the case
    # of non-uniform sampling, but this goes beyond the scope of the example.
    sr = len(list_of_list_of_vars[0][0])/(list_of_list_of_vars[0][0][-1]-list_of_list_of_vars[0][0][0])

    # sampling interval
    ts = 1.0 / sr
    t = np.arange(0, 40, ts)


    X = fft(list_var_interet[1]-np.mean(list_var_interet[1]))
    N = len(X)
    #print(list_var_interet[1])
    n = np.arange(N)
    T = N/sr
    freq = n/T

    # Plot FFT
    fig, ax1 = plt.subplots(figsize=(12, 6))

    ax1.stem(freq, (np.sqrt(np.abs(X) * np.abs(X))) * (np.sqrt(np.abs(X) * np.abs(X)) > 100),
            linefmt='b-', markerfmt=" ", basefmt="-b", label="FFT")
    ax1.set_xlabel('Frequency (Hz)')
    ax1.set_ylabel('FFT Amplitude |X(freq)|', color='b')
    ax1.set_xlim(0.5, 1.5)
    ax1.set_title("FFT of Lift Coefficient")

    # Plot the reference data
    ax2 = ax1.twinx()
    ax2.plot(x_ref, y_ref, 'r', label="Reference Data")
    ax2.set_ylabel('Reference Frequency Spectrum', color='r')

    # Add legends
    fig.legend(loc="upper right", bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)

    plt.show()


spectralAnalysis()



