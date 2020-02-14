import numpy as np
import os
import subprocess
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


################################################################################
# 5 a
################################################################################

print("Running MC simulation for Exercise 5 a. This may take some time ...")
subprocess.call(["./exercise_5_a"])

m = 200
n = 600
l = (m - 1) * (n - 1)

data = np.fromfile("results_5_a.bin", dtype = np.float32)
results_a = data[:l].reshape(m - 1, n - 1, order = "F")
results_s = data[l:2 * l].reshape(m - 1, n - 1, order = "F")
results_b = data[2 * l: 2 * l + 6]
results_sf = data[2 * l + 6:]

try:
    path = os.path.dirname(os.path.abspath(__file__))
except:
    path = "."
plt.style.use(os.path.join(path, "..", "misc", "matplotlib_style.rc"))

f = plt.figure(figsize = (8, 4))
gs = GridSpec(1, 2, width_ratios = [1.0, 0.03])

ax = plt.subplot(gs[0, 0])
xx = np.linspace(0, 60, 600)
yy = np.linspace(0, 20, 200)
img = ax.pcolormesh(xx, yy, results_a / 1e6)

ax.set_ylim([0, 10])
ax.set_xlim([0, 60])
ax.set_xlabel("$x$ [km]")
ax.set_ylabel("$z$ [km]")

ax = plt.subplot(gs[0, 1])
plt.colorbar(img, cax = ax, label = r"Absorbed intensity [$I_{\nu,0}$]")

plt.tight_layout()
f.savefig("results_5_a.png")

################################################################################
# 5 b
################################################################################

print("Running MC simulation for Exercise 5 b. This may take some time ...")
subprocess.call(["./exercise_5_b"])

m = 200
n = 600
l = (m - 1) * (n - 1)

data = np.fromfile("results_5_b.bin", dtype = np.float32)
results_a = data[:l].reshape(m - 1, n - 1, order = "F")
results_s = data[l:2 * l].reshape(m - 1, n - 1, order = "F")
results_b = data[2 * l: 2 * l + 6]
results_sf = data[2 * l + 6:]

f = plt.figure(figsize = (8, 4))
gs = GridSpec(1, 2, width_ratios = [1.0, 0.03])

ax = plt.subplot(gs[0, 0])
xx = np.linspace(0, 60, 600)
yy = np.linspace(0, 20, 200)
img = ax.pcolormesh(xx, yy, results_a / 1e6)

ax.set_ylim([0, 10])
ax.set_xlim([0, 60])
ax.set_xlabel("$x$ [km]")
ax.set_ylabel("$z$ [km]")

ax = plt.subplot(gs[0, 1])
plt.colorbar(img, cax = ax, label = r"Absorbed intensity [$I_{\nu,0}$]")

plt.tight_layout()
f.savefig("results_5_b.png")
