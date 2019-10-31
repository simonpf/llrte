import numpy as np
import os
import subprocess
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

#
# 1a
#

# Run MC
print("Running Monte Carlo simulation for Exercise 2a.")
subprocess.call(["./exercise_2_a"])

# Load matplotlib style
print("Plotting results.")
path = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(path, "..", "misc", "matplotlib_style.rc"))

results_a = []
results_s = []
results_b = []
results_sf = []
for i in range(3):
    n = int(10 ** (i+ 1))
    data = np.fromfile("results_2_a_{}.bin".format(i + 1), dtype = np.float32)
    results_a += [data[:n]]
    results_s += [data[n:2 * n]]
    results_b += [data[2 * n: 2 * n + 6]]
    results_sf += [data[2 * n + 6:]]

#
# Absorbed intensity
#

n_photons = 10000
for i in range(3):
    n = int(10 ** (i + 1))
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    plt.plot(x, results_a[i] / dx, c= "C{}".format(i))
    #plt.plot(x, results_s[i] / dx, c= "C{}".format(i), ls = "--")


#
# Number of scattering events
#

plt.figure()
n_photons = 10000
x = np.arange(11)
for i in range(3):
    plt.bar(x - 0.25 + i / 4.0, results_sf[i] / n_photons, color= "C{}".format(i), width = 1.0 / 4.0)
    #plt.plot(x, results_s[i] / dx, c= "C{}".format(i), ls = "--")


plt.tight_layout()
plt.savefig("results_2_a.png", bbox_inches = "tight", dpi = 300)
