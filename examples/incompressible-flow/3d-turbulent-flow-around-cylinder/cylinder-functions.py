import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.fft import rfft, rfftfreq
from scipy.interpolate import interp1d

def compute_cp_average(
    folder: str,
    theta_array: np.ndarray,
    z_array: np.ndarray,
    cylinder_center: tuple = (8, 25),
    rho: float = 1.0,
    U_inf: float = 1.0,
):
    """
    Compute the time-averaged pressure coefficient (Cp) around a cylinder.

    Args:
        folder (str): Folder containing .vtu files.
        theta_array (np.ndarray): Array of angular positions theta (in radians).
        z_array (np.ndarray): Array of z positions.
        cylinder_center (tuple): (x, y) coordinates of the cylinder center.
        rho (float): Fluid density (default: 1.0).
        U_inf (float): Free-stream velocity (default: 1.0).

    Returns:
        tuple:
            - cp_mean (np.ndarray): Time-averaged Cp as a function of theta.
            - theta_array (np.ndarray): Same input theta array (for convenience).
            - p_inf (float): Reference pressure used for Cp computation.
    """
    np.bool = bool  # Patch for VTK compatibility

    import glob, os

    field="average_pressure"

    vtu_files = sorted(glob.glob(os.path.join(folder, "*.vtu")))
    n_angles = len(theta_array)
    n_timesteps = len(vtu_files)

    r = 0.501
    xc, yc = cylinder_center

    # Determine reference pressure p_inf from the wake
    sim0 = pv.read(vtu_files[-1])
    x_inf = -7.999999
    y_inf = yc
    z_inf = np.mean(z_array)
    p_inf = sim0.sample_over_line([x_inf, y_inf, z_inf], [x_inf, y_inf, z_inf], resolution=1)[field][0]

    print(f"Reference pressure p_inf = {p_inf}")
    # Pressure extraction matrix: rows = time steps, columns = theta angles
    pressure_on_cylinder = np.zeros(n_angles)


    vtu = vtu_files[-1]
    print(f"[{n_timesteps}/{n_timesteps}] Reading: {vtu}")
    sim = pv.read(vtu)
    for j, theta in enumerate(theta_array):
        x = xc + r * np.cos(theta)
        y = yc + r * np.sin(theta)
        p_sum = 0
        for z in z_array:
            point = [x, y, z]
            result = sim.sample_over_line(point, point, resolution=1)
            p_sum += result[field][0]
        pressure_on_cylinder[j] = p_sum / len(z_array)

    # Use only the last N time steps for averaging
    
    cp_mean = (pressure_on_cylinder - p_inf) / (0.5 * rho * U_inf**2)

    return cp_mean, theta_array, p_inf

def compute_strouhal_from_lift(file_path: str, D: float, U: float, method: str = "peak", min_peak_distance: int = 10, purcentage: float = 50):
    """
    Compute the Strouhal number from lift coefficient data.

    Args:
        file_path (str): Path to the force data file.
        D (float): Characteristic length (e.g., cylinder diameter).
        U (float): Free-stream velocity.
        method (str): Method to use: "peak" or "fft".
        min_peak_distance (int): Minimum peak distance for peak method.
        purcentage (float): the purcentage of last time steps to include in the computing in (%)
    Returns:
        float: Estimated Strouhal number.
    """
    # Load data and extract second half
    data = np.genfromtxt(file_path, skip_header=1)
    fraction = int(len(data[:, 0]) / (1/ (purcentage/100)))
    t = data[fraction:, 0]
    cl = data[fraction:, 2]
    cl -= np.mean(cl)

    if method == "peak":
        def smooth(signal, window_size=10):
            return np.convolve(signal, np.ones(window_size)/window_size, mode='same')

        cl_smooth = smooth(cl)
        peaks, _ = find_peaks(cl_smooth, height=0, distance=min_peak_distance)
        t_peaks = t[peaks]
        periods = np.diff(t_peaks)
        T = np.mean(periods)
        f = 1 / T

        plt.plot(t, cl, label="Raw Cl", alpha=0.4)
        plt.plot(t, cl_smooth, label="Smoothed Cl", linewidth=2)
        plt.plot(t_peaks, cl_smooth[peaks], "ro", label="Detected Peaks")
        plt.xlabel("Time (s)")
        plt.ylabel("Cl")
        plt.title("Lift coefficient peak detection for Strouhal number")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    elif method == "fft":
        f_interp = interp1d(t, cl, kind='linear')
        n_steps = len(t) * 20
        x_new = np.linspace(t.min(), t.max(), n_steps)
        y_new = f_interp(x_new)
        dt = (t.max() - t.min()) / n_steps

        Y = rfft(y_new)
        frequencies = rfftfreq(len(y_new), dt)
        amplitudes = np.abs(Y)
        amplitudes[0] = 0  # Ignore DC component
        freq_dominante = frequencies[np.argmax(amplitudes)]
        f = freq_dominante

        plt.plot(t, cl, 'b.', label='Original data')
        plt.plot(x_new, y_new, 'r-', label='Cubic interpolation')
        plt.legend()
        plt.title("Cubic interpolation of Cl")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        plt.loglog(frequencies, amplitudes)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Amplitude")
        plt.title("FFT spectrum of interpolated Cl")
    
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    else:
        raise ValueError("Method must be either 'peak' or 'fft'")

    St = f * D / U
    return St



def compute_drag_coefficient(file_path: str, D: float, L: float, rho: float = 1.0, U: float = 1.0, n_final_timesteps: int = 2000):
    """
    Compute the mean drag coefficient Cd from a .dat force file.

    Args:
        file_path (str): Path to the force data file.
        D (float): Diameter of the cylinder.
        L (float): Length of the cylinder.
        rho (float): Fluid density.
        U (float): Free-stream velocity.
        n_final_timesteps (int): Number of last time steps to consider.

    Returns:
        float: Mean drag coefficient.
    """
    data = np.genfromtxt(file_path, skip_header=1)
    t = data[-n_final_timesteps:, 0]
    drag = data[-n_final_timesteps:, 1]
    cd = drag / (0.5 * rho * U**2 * D * L)
    cd_mean = np.mean(cd)



    return cd_mean
