import numpy as np
import os
import subprocess
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

#
# 2 a
#

# Run MC
print("Running Monte Carlo simulation for Exercise 2a ...")
subprocess.call(["./exercise_2_a"])

# Load matplotlib style
print("Plotting results ...\n")
path = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(path, "..", "misc", "matplotlib_style.rc"))

################################################################################
# Grid resolution
################################################################################

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

f = plt.figure(figsize = (20, 5))
gs = GridSpec(1, 4, width_ratios = [1.0, 1.0, 1.0, 1.0])
handles = []

#
# Absorbed intensity
#

ax = plt.subplot(gs[0])
n_photons = 10000
for i in range(3):
    n = int(10 ** (i + 1))
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    handles += ax.plot(x, results_a[i] / dx / n_photons, c= "C{}".format(i))
    #plt.plot(x, results_s[i] / dx, c= "C{}".format(i), ls = "--")
ax.set_title("(a) Absorbed intensity", loc = "left")

ax.set_ylim([0, 3e-5])
ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Absorbed intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")

#
# Scattered intensity
#

ax = plt.subplot(gs[1])
n_photons = 10000
for i in range(2, -1, -1):
    n = int(10 ** (i + 1))
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    ax.scatter(x, results_s[i] / n_photons / dx, c= "C{}".format(i))

ax.set_ylim([0, 1.5e-4])
ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Scattered intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")
ax.set_title("(b) Scattered intensity", loc = "left")


#
# Number of scattering events
#

ax = plt.subplot(gs[2])
x = np.arange(11)
for i in range(3):
    plt.bar(x - 0.25 + i / 4.0, results_sf[i] / n_photons, color= "C{}".format(i), width = 1.0 / 4.0)
    #plt.plot(x, results_s[i] / dx, c= "C{}".format(i), ls = "--")
ax.set_xticks(np.arange(11))
ax.set_xticklabels([str(i) for i in range(11)])
ax.set_xlabel("Number of times scattered")
ax.set_ylabel("Frequency")
ax.set_title("(c) Number of scattering events", loc = "left")

plt.tight_layout()

ax_l = plt.subplot(gs[-1])
labels = [r"$n_\text{{grid}}$ = {}".format(n) for n in [1000, 100, 10][::-1]]
ax_l.set_axis_off()
ax_l.legend(handles = handles, labels = labels, loc = "center")

f.savefig("results_2_a_1.png", bbox_inches = "tight", dpi = 300)

print("Absorbed intensity for varying grid size:")
print("\t Left boundary:  {:.4} / {:.4} / {:.4}".format(*[r[0] / n_photons for r in results_b]))
print("\t Domain:         {:.4} / {:.4} / {:.4}".format(*[r.sum() / n_photons for r in results_a]))
print("\t Right boundary: {:.4} / {:.4} / {:.4}\n\n".format(*[r[1] / n_photons for r in results_b]))

################################################################################
# Photon counts
################################################################################

results_a = []
results_s = []
results_b = []
results_sf = []
for i in range(3):
    n = 100
    data = np.fromfile("results_2_a_{}.bin".format(i + 4), dtype = np.float32)
    results_a += [data[:n]]
    results_s += [data[n:2 * n]]
    results_b += [data[2 * n: 2 * n + 6]]
    results_sf += [data[2 * n + 6:]]

f = plt.figure(figsize = (20, 5))
gs = GridSpec(1, 4, width_ratios = [1.0, 1.0, 1.0, 1.0])
handles = []

#
# Absorbed intensity
#

ax = plt.subplot(gs[0])
n = 100
for i in range(3):
    n_photons = int(10 ** ((i + 1) * 2))
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    handles += ax.plot(x, results_a[i] / dx / n_photons, c= "C{}".format(i))
ax.set_title("(a) Absorbed intensity", loc = "left")
ax.set_ylim([0, 3e-5])

ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Absorbed intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")

#
# Scattered intensity
#

ax = plt.subplot(gs[1])
n = 100
for i in range(2, -1, -1):
    n_photons = int(10 ** ((i + 1) * 2))
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    ax.scatter(x, results_s[i] / dx / n_photons, c= "C{}".format(i))

ax.set_ylim([0, 1.5e-4])
ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Scattered intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")
ax.set_title("(b) Scattered intensity", loc = "left")


#
# Number of scattering events
#

ax = plt.subplot(gs[2])
x = np.arange(11)
for i in range(3):
    n_photons = int(10 ** ((i + 1) * 2))
    plt.bar(x - 0.25 + i / 4.0, results_sf[i] / n_photons, color= "C{}".format(i), width = 1.0 / 4.0)
    #plt.plot(x, results_s[i] / dx, c= "C{}".format(i), ls = "--")
ax.set_xticks(np.arange(11))
ax.set_xticklabels([str(i) for i in range(11)])
ax.set_xlabel("Number of times scattered")
ax.set_ylabel("Frequency")
ax.set_title("(c) Number of scattering events", loc = "left")


ax_l = plt.subplot(gs[-1])
labels = [r"$n_\text{{photons}}$ = {}".format(n) for n in [10, 100, 1000]]
ax_l.set_axis_off()
ax_l.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()
f.savefig("results_2_a_2.png", bbox_inches = "tight", dpi = 300)

print("Absorbed intensity for varying photon counts:\n")
print("\t Left boundary:  {:.4} / {:.4} / {:.4}".format(*[r[0] / (10 ** ((i + 1) * 2)) for i,r in enumerate(results_b)]))
print("\t Domain:         {:.4} / {:.4} / {:.4}".format(*[r.sum() / (10 ** ((i + 1) * 2)) for i,r in enumerate(results_a)]))
print("\t Right boundary: {:.4} / {:.4} / {:.4}\n\n".format(*[r[1] / (10 ** ((i + 1) * 2)) for i,r in enumerate(results_b)]))

################################################################################
# Optical depth
################################################################################

results_a = []
results_s = []
results_b = []
results_sf = []
for i in range(3):
    n = 100
    data = np.fromfile("results_2_a_{}.bin".format(i + 7), dtype = np.float32)
    results_a += [data[:n]]
    results_s += [data[n:2 * n]]
    results_b += [data[2 * n: 2 * n + 6]]
    results_sf += [data[2 * n + 6:]]

f = plt.figure(figsize = (20, 5))
gs = GridSpec(1, 4, width_ratios = [1.0, 1.0, 1.0, 1.0])
handles = []

#
# Absorbed intensity
#

ax = plt.subplot(gs[0])
n = 100
n_photons = 10000
for i in range(3):
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    handles += ax.plot(x, results_a[i] / dx / n_photons, c= "C{}".format(i))
ax.set_title("(a) Absorbed intensity", loc = "left")
ax.set_ylim([0, 1e-4])

ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Absorbed intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")

#
# Scattered intensity
#

ax = plt.subplot(gs[1])
n = 100
n_photons = 10000
for i in range(2, -1, -1):
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    ax.scatter(x, results_s[i] / dx / n_photons, c= "C{}".format(i))

ax.set_ylim([0, 5e-4])
ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Scattered intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")
ax.set_title("(b) Scattered intensity", loc = "left")


#
# Number of scattering events
#

ax = plt.subplot(gs[2])
x = np.arange(11)
n = 100
for i in range(3):
    plt.bar(x - 0.25 + i / 4.0, results_sf[i] / n_photons, color= "C{}".format(i), width = 1.0 / 4.0)
    #plt.plot(x, results_s[i] / dx, c= "C{}".format(i), ls = "--")
ax.set_xticks(np.arange(11))
ax.set_xticklabels([str(i) for i in range(11)])
ax.set_xlabel("Number of times scattered")
ax.set_ylabel("Frequency")
ax.set_title("(c) Number of scattering events", loc = "left")


ax_l = plt.subplot(gs[-1])
labels = [r"$\tau = {}$".format(n) for n in [1, 2, 3]]
ax_l.set_axis_off()
ax_l.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()
f.savefig("results_2_a_3.png", bbox_inches = "tight", dpi = 300)

n = 100
print("Absorbed intensity for varying optical depth:\n")

print("\t Left boundary:  {:.4} / {:.4} / {:.4}".format(*[r[0] / n_photons for i,r in enumerate(results_b)]))
print("\t Domain:         {:.4} / {:.4} / {:.4}".format(*[r.sum() / n_photons for i,r in enumerate(results_a)]))
print("\t Right boundary: {:.4} / {:.4} / {:.4}\n\n".format(*[r[1] / n_photons for i,r in enumerate(results_b)]))

################################################################################
# Forward-to-backscattering ratio
################################################################################

results_a = []
results_s = []
results_b = []
results_sf = []
for i in range(3):
    n = 100
    data = np.fromfile("results_2_a_{}.bin".format(i + 10), dtype = np.float32)
    results_a += [data[:n]]
    results_s += [data[n:2 * n]]
    results_b += [data[2 * n: 2 * n + 6]]
    results_sf += [data[2 * n + 6:]]

f = plt.figure(figsize = (20, 5))
gs = GridSpec(1, 4, width_ratios = [1.0, 1.0, 1.0, 1.0])
handles = []

#
# Absorbed intensity
#

ax = plt.subplot(gs[0])
n = 100
n_photons = 10000
for i in range(3):
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    handles += ax.plot(x, results_a[i] / dx / n_photons, c= "C{}".format(i))
ax.set_title("(a) Absorbed intensity", loc = "left")
ax.set_ylim([0, 3e-5])

ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Absorbed intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")

#
# Scattered intensity
#

ax = plt.subplot(gs[1])
n = 100
for i in range(2, -1, -1):
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    ax.scatter(x, results_s[i] / dx / n_photons, c= "C{}".format(i))

ax.set_ylim([0, 1.5e-4])
ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Scattered intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")
ax.set_title("(b) Scattered intensity", loc = "left")


#
# Number of scattering events
#

ax = plt.subplot(gs[2])
x = np.arange(11)
n = 100
for i in range(3):
    plt.bar(x - 0.25 + i / 4.0, results_sf[i] / n_photons, color= "C{}".format(i), width = 1.0 / 4.0)
    #plt.plot(x, results_s[i] / dx, c= "C{}".format(i), ls = "--")
ax.set_xticks(np.arange(11))
ax.set_xticklabels([str(i) for i in range(11)])
ax.set_xlabel("Number of times scattered")
ax.set_ylabel("Frequency")
ax.set_title("(c) Number of scattering events", loc = "left")


ax_l = plt.subplot(gs[-1])
labels = [r"$VR = {}$".format(n) for n in [0.2, 0.5, 0.8]]
ax_l.set_axis_off()
ax_l.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()
f.savefig("results_2_a_4.png", bbox_inches = "tight", dpi = 300)

n = 100
print("Absorbed intensity for varying forward-to-backscattering ratio:\n")

print("\t Left boundary:  {:.4} / {:.4} / {:.4}".format(*[r[0] / n_photons for i,r in enumerate(results_b)]))
print("\t Domain:         {:.4} / {:.4} / {:.4}".format(*[r.sum() / n_photons for i,r in enumerate(results_a)]))
print("\t Right boundary: {:.4} / {:.4} / {:.4}\n\n".format(*[r[1] / n_photons for i,r in enumerate(results_b)]))

################################################################################
# Single scattering albedo
################################################################################

results_a = []
results_s = []
results_b = []
results_sf = []
for i in range(3):
    n = 100
    data = np.fromfile("results_2_a_{}.bin".format(i + 13), dtype = np.float32)
    results_a += [data[:n]]
    results_s += [data[n:2 * n]]
    results_b += [data[2 * n: 2 * n + 6]]
    results_sf += [data[2 * n + 6:]]

f = plt.figure(figsize = (20, 5))
gs = GridSpec(1, 4, width_ratios = [1.0, 1.0, 1.0, 1.0])
handles = []

#
# Absorbed intensity
#

ax = plt.subplot(gs[0])
n = 100
n_photons = 10000
for i in range(3):
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    handles += ax.plot(x, results_a[i] / dx / n_photons, c= "C{}".format(i))
ax.set_title("(a) Absorbed intensity", loc = "left")
ax.set_ylim([0, 1e-4])

ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Absorbed intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")

#
# Scattered intensity
#

ax = plt.subplot(gs[1])
n = 100
for i in range(2, -1, -1):
    x = np.linspace(0, 10, n + 1)
    x = x[:-1] + 0.5 * np.diff(x)[0]
    dx = 10e3 / n
    ax.scatter(x, results_s[i] / dx / n_photons, c= "C{}".format(i))

ax.set_ylim([0, 1.5e-4])
ax.set_xlabel("Distance [km]")
ax.set_ylabel(r"Scattered intensity [$I_{\nu,0}\ \unit{km}^{-1}$]")
ax.set_title("(b) Scattered intensity", loc = "left")


#
# Number of scattering events
#

ax = plt.subplot(gs[2])
x = np.arange(11)
n = 100
for i in range(3):
    plt.bar(x - 0.25 + i / 4.0, results_sf[i] / n_photons, color= "C{}".format(i), width = 1.0 / 4.0)
    #plt.plot(x, results_s[i] / dx, c= "C{}".format(i), ls = "--")
ax.set_xticks(np.arange(11))
ax.set_xticklabels([str(i) for i in range(11)])
ax.set_xlabel("Number of times scattered")
ax.set_ylabel("Frequency")
ax.set_title("(c) Number of scattering events", loc = "left")


ax_l = plt.subplot(gs[-1])
labels = [r"$SSA = {}$".format(n) for n in [0.2, 0.5, 0.8]]
ax_l.set_axis_off()
ax_l.legend(handles = handles, labels = labels, loc = "center")

plt.tight_layout()
f.savefig("results_2_a_5.png", bbox_inches = "tight", dpi = 300)

n = 100
print("Absorbed intensity for varying single-scattering albedo:\n")

print("\t Left boundary:  {:.4} / {:.4} / {:.4}".format(*[r[0] / n_photons for i,r in enumerate(results_b)]))
print("\t Domain:         {:.4} / {:.4} / {:.4}".format(*[r.sum() / n_photons for i,r in enumerate(results_a)]))
print("\t Right boundary: {:.4} / {:.4} / {:.4}".format(*[r[1] / n_photons for i,r in enumerate(results_b)]))
