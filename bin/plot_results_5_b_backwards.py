import numpy as np
import os
import subprocess
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xarray

################################################################################
# 5 a
################################################################################

data = xarray.open_dataset("exercise_5_b_backwards.nc")
y = data.data[:, 0, 0, 0]
x = np.linspace(0, 60, y.size)
plt.plot(x, y)
plt.show()
