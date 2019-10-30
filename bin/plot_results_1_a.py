import numpy as np
import os
import subprocess
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

#
# 1a
#

# Run MC
print("Running Monte Carlo simulation for Exercise 1a.")
subprocess.call(["./exercise_1_a"])

# Load matplotlib style
print("Plotting results.")
path = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(path, "..", "misc", "matplotlib_style.rc"))

results = []
for i in range(3):
    results += [np.fromfile("results_1_a_{}.bin".format(i + 1), dtype = np.float32)]

# Results of Beer-Lambert law
x = np.linspace(0, 10e3, 1001)
sigma = np.zeros(x.shape)
sigma[:] = 1e-3
tau = cumtrapz(sigma, x = x)
plt.plot(x[1:] / 1e3, np.exp(-tau), ls = "--", c = "k", label = "Beer-Lambert Law")

# Plot data
n_photons = 10000
for d, n in zip(results[::-1], [10, 100, 1000][::-1]):
    x = np.linspace(0, 10e3, n + 1)[:-1]
    plt.scatter(x / 1e3,
                d / n_photons,
                label = r"$n_\text{{grid}}$ = {}".format(n),
                s = 5,
                marker = "x")
plt.legend()

plt.ylim([0, 1])
plt.xlim([0, 10])
plt.ylabel(r"$I_\nu$ [$I_{\nu, 0}$]")
plt.xlabel("Distance [km]")
plt.tight_layout()
filename = "results_1_a_1.png"
print("Saving results to {}".format(filename))
plt.savefig(filename, bbox_inches = "tight", dpi = 300)

#
# 1a
#

plt.figure()
results = []
for i in range(3):
    results += [np.fromfile("results_1_a_{}.bin".format(i + 4), dtype = np.float32)]

# Results of Beer-Lambert law
x = np.linspace(0, 10e3, 1001)
sigma = np.zeros(x.shape)
sigma[:] = 1e-3
tau = cumtrapz(sigma, x = x)
plt.plot(x[1:] / 1e3, np.exp(-tau), ls = "--", c = "k", label = "Beer-Lambert Law")

# Plot data
for d, n in zip(results, [100, 10000, 1000000]):
    x = np.linspace(0, 10e3, 101)[:-1]
    plt.scatter(x / 1e3,
                d / n,
                label = r"$n_\text{{photons}}$ = {}".format(n),
                s = 5,
                marker = "x")
plt.legend()

plt.ylim([0, 1])
plt.xlim([0, 10])
plt.ylabel(r"$I_\nu$ [$I_{\nu, 0}$]")
plt.xlabel("Distance [km]")
plt.tight_layout()
filename = "results_1_a_2.png"
print("Saving results to {}".format(filename))
plt.savefig(filename, bbox_inches = "tight", dpi = 300)

plt.figure()
