import numpy as np
import os
import subprocess
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


################################################################################
# 4 b 1
################################################################################

print("Running MC simulation for Exercise 4 b. This may take some time ...")
proc = subprocess.Popen("./exercise_4_b", stdout = subprocess.PIPE)
results = np.loadtxt(proc.stdout)

try:
    path = os.path.dirname(os.path.abspath(__file__))
except:
    path = "."
plt.style.use(os.path.join(path, "..", "misc", "matplotlib_style.rc"))

f = plt.figure(figsize = (8, 4))
gs = GridSpec(1, 2, width_ratios = [1.0, 0.2])

g = results[:, 0]
handles = []
labels = []

ax = plt.subplot(gs[0, 0])
handles += ax.plot(g, results[:, 1], marker = "o")
labels += ["Upper boundary"]
handles += ax.plot(g, results[:, 2], marker = "o")
labels += ["Medium"]
handles += ax.plot(g, results[:, 3], marker = "o")
labels += ["Lower surface"]
handles += ax.plot(g, results[:, 1:].sum(axis = -1), c = "grey", ls = "--")
labels += ["Total"]
ax.set_ylabel(r"Absorbed intensity [$I_{\nu, 0}$]")
ax.set_xlabel("Asymmetry parameter")
ax.set_xlim(-2.0, 1.0)
ax.set_ylim(0.0, 1.1)

ax = plt.subplot(gs[0, 1])
ax.set_axis_off()
ax.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()
f.savefig("results_4_b.png")


################################################################################
# Henyey-Greenstein phase function
################################################################################

f = plt.figure(figsize = (8, 4))
gs = GridSpec(1, 2, width_ratios = [1.0, 0.2])

g = results[:, 0]
handles = []
labels = []

def hg(theta, g):
    theta = np.pi * theta / 180.0
    return (1 - g ** 2) / (1 + g ** 2 - 2 * g * np.cos(theta)) ** 1.5

theta = np.linspace(0, 180, 101)
g = np.array([-0.8, -0.4, 0.0, 0.4, 0.8])
ps = hg(theta.reshape(1, -1), g.reshape(-1, 1))
ps /= np.trapz(ps, theta, axis = -1).reshape(-1, 1) * 4 * np.pi

from matplotlib.cm import magma
cs = [magma(i / 5) for i in range(ps.shape[0])]

ax = plt.subplot(gs[0, 0])

for i in range(ps.shape[0]):
    handles += ax.plot(theta, ps[i, :], c = cs[i])
    labels += ["g = {}".format(g[i])]
ax.set_xlabel(r"$\theta$ [$^\circ$]")
ax.set_ylabel(r"$p(\theta)$")
ax.set_yscale("log")
ax.set_xlim(0, 180)

ax = plt.subplot(gs[0, 1])
ax.set_axis_off()
ax.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()
f.savefig("henyey_greenstein.png")
