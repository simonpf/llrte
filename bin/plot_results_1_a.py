import numpy as np
import os
import subprocess
import scipy as sp
import xarray
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

################################################################################
# 1a
################################################################################

# Run MC
print("Running Monte Carlo simulation for Exercise 1 a.")
subprocess.call(["./exercise_1_a"])

# Load matplotlib style
print("Plotting results.")
path = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(path, "..", "misc", "matplotlib_style.rc"))

results = []
for i in range(3):
    data = xarray.open_dataset("results_1_a_{}.nc".format(i + 1))
    results += [np.array(data.irradiance)]


# Results of Beer-Lambert law
x = np.linspace(0, 10e3, 1001)
sigma = np.zeros(x.shape)
sigma[:] = 1e-4
tau = cumtrapz(sigma, x = x)

f = plt.figure(figsize = (10, 5))
gs = GridSpec(1, 2, width_ratios = [1.0, 1.0])
ax = plt.subplot(gs[0])
handles = []
handles += ax.plot(x[1:] / 1e3, np.exp(-tau), ls = "--", c = "k", label = "Beer-Lambert Law")

# Plot data
n_photons = 10000
for d in results[::-1]:
    n = d.shape[0]
    x = np.linspace(0, 10e3, n - 1)[:-2]
    handles += [ax.scatter(x / 1e3,
                           d[2:-1, 1, 1] / n_photons,
                           label = r"$n_\text{{grid}}$ = {}".format(n),
                           marker = "x")]

ax_l = plt.subplot(gs[1])
labels = ["Absorption law"] + [r"$n_\text{{grid}}$ = {}".format(n) for n in [1000, 100, 10]]
ax_l.set_axis_off()
ax_l.legend(handles = handles, labels = labels, loc = "center")

ax.set_ylim([0, 1])
ax.set_xlim([0, 10])
ax.set_ylabel(r"$I_\nu$ [$I_{\nu, 0}$]")
ax.set_xlabel("Distance [km]")
plt.tight_layout()
filename = "results_1_a_1.png"
print("Saving results to {}".format(filename))
f.savefig(filename, bbox_inches = "tight", dpi = 300)

################################################################################
# 1b
################################################################################

results = []
for i in range(3):
    data = xarray.open_dataset("results_1_a_{}.nc".format(i + 4))
    results += [np.array(data.irradiance)]

# Results of Beer-Lambert law
x = np.linspace(0, 10e3, 1001)
sigma = np.zeros(x.shape)
sigma[:] = 1e-4
tau = cumtrapz(sigma, x = x)

f = plt.figure(figsize = (10, 5))
gs = GridSpec(1, 2, width_ratios = [1.0, 1.0])
ax = plt.subplot(gs[0])
handles = []
handles += ax.plot(x[1:] / 1e3, np.exp(-tau), ls = "--", c = "k", label = "Beer-Lambert Law")

# Plot data
for d, n in zip(results, [100, 10000, 1000000]):
    x = np.linspace(0, 10e3, 101)[:-2]
    handles += [ax.scatter(x / 1e3,
                           d[2:-1, 1, 1] / n,
                           marker = "x",
                           label = r"$n_\text{{photons}}$ = {}".format(n))]

ax_l = plt.subplot(gs[1])
labels = ["Absorption law"] + [r"$n_\text{{photons}}$ = {}".format(n) for n in [100, 10000, 1000000]]
ax_l.set_axis_off()
ax_l.legend(handles = handles, labels = labels, loc = "center")

ax.set_ylim([0, 1])
ax.set_xlim([0, 10])
ax.set_ylabel(r"$I_\nu$ [$I_{\nu, 0}$]")
ax.set_xlabel("Distance [km]")
plt.tight_layout()
filename = "results_1_a_2.png"
print("Saving results to {}".format(filename))
f.savefig(filename, bbox_inches = "tight", dpi = 300)
