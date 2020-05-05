import numpy as np
import matplotlib.pyplot as plt
import xarray

#
# Lambertian surface
#

data = xarray.open_dataset("surface_lambertian.nc")
out = np.array(data.outgoing_directions)
outp = np.array(data.outgoing_positions)
inc = np.array(data.incoming_directions)
incp = np.array(data.incoming_positions)

thetas = np.arccos(-out[:, 0])
phis = np.arctan2(out[:, 1], out[:, 2])

f, axs = plt.subplots(1, 2)
ax = axs[0]
ax.hist(thetas, normed=True)
ax.set_ylabel("Frequency")
ax.set_xlabel(r"$\theta$")
ax = axs[1]
ax.hist(phis, normed=True)
ax.set_xlabel(r"$\phi$")
plt.tight_layout()
plt.savefig("lambertian_angles.png")

#
# Specular surface
#

data = xarray.open_dataset("surface_specular.nc")
out = np.array(data.outgoing_directions)
outp = np.array(data.outgoing_positions)
inc = np.array(data.incoming_directions)
incp = np.array(data.incoming_positions)

thetas = np.arccos(-out[:, 0])
phis = np.arctan2(out[:, 1], out[:, 2])

f, axs = plt.subplots(1, 2)
ax = axs[0]
ax.hist(thetas, normed=True)
ax.set_ylabel("Frequency")
ax.set_xlabel(r"$\theta$")
ax = axs[1]
ax.hist(phis, normed=True)
ax.set_xlabel(r"$\phi$")
plt.tight_layout()
plt.savefig("specular_angles.png")
