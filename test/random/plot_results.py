import numpy as np
import matplotlib.pyplot as plt
import xarray

data = xarray.open_dataset("file.nc")

data_full = np.array(data.data_full)
bins = np.linspace(0.0, np.pi, 101)
y, x = np.histogram(data_full.ravel(), bins=bins, density=True)
x = 0.5 * (x[1:] + x[:-1])
plt.plot(x, y, zorder=1)
y = np.sin(x) * 0.5
plt.plot(x, y, c="k", ls="--", zorder=0)
plt.savefig("results_full.pdf")

plt.clf()
data_half = np.array(data.data_half)
bins = np.linspace(0.0, np.pi, 101)
y, x = np.histogram(data_half.ravel(), bins=bins, density=True)
x = 0.5 * (x[1:] + x[:-1])
plt.plot(x, y, zorder=1)
y = np.sin(x)
plt.plot(x, y, c="k", ls="--", zorder=0)
plt.savefig("results_half.pdf")
