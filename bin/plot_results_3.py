import numpy as np
import os
import subprocess
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


################################################################################
# 3 b 1
################################################################################

print("Running MC simulation for Exercise 3 b.1. This may take some time ...")
proc = subprocess.Popen("./exercise_3_b_1", stdout = subprocess.PIPE)
results = np.loadtxt(proc.stdout)

try:
    path = os.path.dirname(os.path.abspath(__file__))
except:
    path = "."
plt.style.use(os.path.join(path, "..", "misc", "matplotlib_style.rc"))

f = plt.figure(figsize = (8, 4))
gs = GridSpec(1, 2, width_ratios = [1.0, 0.2])

albedo = np.linspace(0.0, 1.0, 11)
handles = []
labels = []

ax = plt.subplot(gs[0, 0])
handles += ax.plot(albedo, results[:, 0], marker = "o")
labels += ["Upper boundary"]
handles += ax.plot(albedo, results[:, 1], marker = "o")
labels += ["Medium"]
handles += ax.plot(albedo, results[:, 2], marker = "o")
labels += ["Lower surface"]
handles += ax.plot(albedo, results.sum(axis = -1), c = "grey", ls = "--")
labels += ["Total"]
ax.set_ylabel(r"Absorbed intensity [$I_{\nu, 0}$]")
ax.set_xlabel("Surface albedo")
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.1)

ax = plt.subplot(gs[0, 1])
ax.set_axis_off()
ax.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()
f.savefig("results_3_b_1.png")

################################################################################
# 3 b 2
################################################################################

print("Running MC simulation for Exercise 3 b.2. This may take some time ...")
proc = subprocess.Popen("./exercise_3_b_2", stdout = subprocess.PIPE)
results = np.loadtxt(proc.stdout)

f = plt.figure(figsize = (8, 4))
gs = GridSpec(1, 2, width_ratios = [1.0, 0.2])

albedo = np.linspace(0.0, 1.0, 11)
handles = []
labels = []

ax = plt.subplot(gs[0, 0])
handles += ax.plot(albedo, results[:, 0], marker = "o")
labels += ["Upper boundary"]
handles += ax.plot(albedo, results[:, 1], marker = "o")
labels += ["Medium"]
handles += ax.plot(albedo, results[:, 2], marker = "o")
labels += ["Lower surface"]
handles += ax.plot(albedo, results.sum(axis = -1), c = "grey", ls = "--")
labels += ["Total"]
ax.set_ylabel(r"Absorbed intensity [$I_{\nu, 0}$]")
ax.set_xlabel("Single scattering albedo")
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.1)

ax = plt.subplot(gs[0, 1])
ax.set_axis_off()
ax.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()
f.savefig("results_3_b_2.png")

################################################################################
# 3 b 3
################################################################################

print("Running MC simulation for Exercise 3 b.3. This may take some time ...")
proc = subprocess.Popen("./exercise_3_b_3", stdout = subprocess.PIPE)
results = np.loadtxt(proc.stdout)

f = plt.figure(figsize = (8, 4))
gs = GridSpec(1, 2, width_ratios = [1.0, 0.2])

ods = np.linspace(0.2, 2.0, 10)
handles = []
labels = []

ax = plt.subplot(gs[0, 0])
handles += ax.plot(ods, results[:, 0], marker = "o")
labels += ["Upper boundary"]
handles += ax.plot(ods, results[:, 1], marker = "o")
labels += ["Medium"]
handles += ax.plot(ods, results[:, 2], marker = "o")
labels += ["Lower surface"]
handles += ax.plot(ods, results.sum(axis = -1), c = "grey", ls = "--")
labels += ["Total"]
ax.set_ylabel(r"Absorbed intensity [$I_{\nu, 0}$]")
ax.set_xlabel("Optical depth")
ax.set_xlim(0.2, 2.0)
ax.set_ylim(0.0, 1.1)

ax = plt.subplot(gs[0, 1])
ax.set_axis_off()
ax.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()
f.savefig("results_3_b_3.png")
