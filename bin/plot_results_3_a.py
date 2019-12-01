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
print("Running Monte Carlo simulation for Exercise 3 a ...")
subprocess.call(["./exercise_3_a"])

# Load matplotlib style
print("Plotting results ...\n")
path = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(path, "..", "misc", "matplotlib_style.rc"))

################################################################################
# Grid resolution
################################################################################

data = np.fromfile("results_3_a.bin", dtype = np.float32)
n = 100 * 100
results_a = data[:n].reshape(100, 100)
results_s = data[n:2 * n].reshape(100, 100)
results_b = data[2 * n : 2 * n + 6]
results_sf = data[2 * n + 6 :]
