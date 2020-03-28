import numpy as np
import os
import subprocess
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xarray

################################################################################
# 5 3d
################################################################################

data = xarray.open_dataset("exercise_5_3d.nc")
data = np.array(data.data[:, :, 0, 0])
x = np.linspace(-10, 10, data.shape[0] + 1)
y = np.linspace(-10, 10, data.shape[0] + 1)
plt.pcolormesh(x, y, data)
