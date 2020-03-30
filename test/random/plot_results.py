import numpy as np
import matplotlib.pyplot as plt
import xarray

#
# Exponential distribution
#

plt.clf()
data = xarray.open_dataset("data_exp.nc")
data = np.array(data.data)
data_cuda = xarray.open_dataset("data_exp_cuda.nc")
data_cuda = np.array(data_cuda.data)

bins = np.linspace(0.0, 10, 101)
y, x = np.histogram(data.ravel(), bins=bins, density=True)
y_cuda, x = np.histogram(data_cuda.ravel(), bins=bins, density=True)
x = 0.5 * (x[1:] + x[:-1])
plt.plot(x, y, label = "Host")
plt.plot(x, y_cuda, label = "GPU")
plt.legend()
plt.savefig("results_exp.pdf")

#
# Phit
#

plt.clf()
data = xarray.open_dataset("data_phi.nc")
data = np.array(data.data)
data_cuda = xarray.open_dataset("data_phi_cuda.nc")
data_cuda = np.array(data_cuda.data)

bins = np.linspace(0.0, 10, 101)
y, x = np.histogram(data.ravel(), bins=bins, density=True)
y_cuda, x = np.histogram(data_cuda.ravel(), bins=bins, density=True)
x = 0.5 * (x[1:] + x[:-1])
plt.plot(x, y, label = "Host")
plt.plot(x, y_cuda, label = "GPU")
plt.legend()
plt.savefig("results_phi.pdf")

#
# Theta
#

plt.clf()
data = xarray.open_dataset("data_theta.nc")
data = np.array(data.data)
data_cuda = xarray.open_dataset("data_theta_cuda.nc")
data_cuda = np.array(data_cuda.data)

bins = np.linspace(0.0, 10, 101)
y, x = np.histogram(data.ravel(), bins=bins, density=True)
y_cuda, x = np.histogram(data_cuda.ravel(), bins=bins, density=True)
x = 0.5 * (x[1:] + x[:-1])
plt.plot(x, y, label = "Host")
plt.plot(x, y_cuda, label = "GPU")
plt.legend()
plt.savefig("results_theta.pdf")
