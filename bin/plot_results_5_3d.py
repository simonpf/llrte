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
x = np.array(data.x)
y = np.array(data.y)
data = np.array(data.data[:, :, 0, 0])
plt.pcolormesh(x, y, data)
plt.show()
