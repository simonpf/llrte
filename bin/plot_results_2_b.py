import numpy as np
import os
import subprocess
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

#
# 2 b
#

# Run MC
print("Running Monte Carlo simulation for Exercise 2b ...")
subprocess.call(["./exercise_2_b"])

# Load matplotlib style
print("Plotting results ...\n")
path = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(path, "..", "misc", "matplotlib_style.rc"))

################################################################################
# Grid resolution
################################################################################

n_photons = 10000
n = 100

data = np.fromfile("results_2_b.bin", dtype = np.float32)
results_a = data[:n]
results_s = data[n:2 * n]
results_b = data[2 * n: 2 * n + 6]
results_sf = data[2 * n + 6:]

data_r = np.fromfile("results_2_a_10.bin", dtype = np.float32)
results_a_r = data_r[:n]
results_s_r = data_r[n:2 * n]
results_b_r = data_r[2 * n: 2 * n + 6]
results_sf_r = data_r[2 * n + 6:]

f = plt.figure(figsize = (20, 5))
gs = GridSpec(1, 4, width_ratios = [1.0, 1.0, 1.0, 1.0])
handles = []

#
# Absorbed intensity
#

ax = plt.subplot(gs[0])
x = np.linspace(0, 10, n + 1)
x = x[:-1] + 0.5 * np.diff(x)[0]
dx = 10e3 / n
handles += ax.plot(x, results_a / dx / n_photons)
handles += ax.plot(x, results_a_r / dx / n_photons)

ax.set_title("(a) Absorbed intensity", loc = "left")
ax.set_ylim([0, 3e-5])
ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Absorbed intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")

#
# Scattered intensity
#

ax = plt.subplot(gs[1])
x = np.linspace(0, 10, n + 1)
x = x[:-1] + 0.5 * np.diff(x)[0]
dx = 10e3 / n
ax.scatter(x, results_s / n_photons / dx)
ax.scatter(x, results_s_r / dx / n_photons)

ax.set_ylim([0, 1.5e-4])
ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Scattered intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")
ax.set_title("(b) Scattered intensity", loc = "left")


#
# Number of scattering events
#

ax = plt.subplot(gs[2])
x = np.arange(11)
plt.bar(x - 0.33, results_sf / n_photons, width = 0.4)
plt.bar(x, results_sf_r / n_photons, width = 0.4)

ax.set_xticks(np.arange(11))
ax.set_xticklabels([str(i) for i in range(11)])
ax.set_xlabel("Number of times scattered")
ax.set_ylabel("Frequency")
ax.set_title("(c) Number of scattering events", loc = "left")


ax_l = plt.subplot(gs[-1])
labels = ["Heterogeneous", "Homogeneous"]
ax_l.set_axis_off()
ax_l.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()


f.savefig("results_2_b.png", bbox_inches = "tight", dpi = 300)
