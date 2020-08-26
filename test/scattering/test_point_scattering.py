import numpy as np
import scipy as sp
import scipy.integrate
import os
import netCDF4
from gridded_optical_properties import PointScattering

import matplotlib.pyplot as plt

class ScatteringData:
    def __init__(self, filename):
        handle = netCDF4.Dataset(filename)
        self.azimuth_angles_incoming = handle["azimuth_angles_incoming"][:].data
        self.zenith_angles_incoming = handle["zenith_angles_incoming"][:].data
        self.azimuth_angles_scattering = handle["azimuth_angles_scattering"][:].data
        self.zenith_angles_scattering = handle["zenith_angles_scattering"][:].data
        self.phase_matrix = handle["phase_matrix"][:].data
        v = handle["absorption_vector"][:].data
        self.absorption_vector = v.reshape((1, 1,) + v.shape)
        v = handle["extinction_matrix"][:].data
        self.extinction_matrix = v.reshape((1, 1,) + v.shape)
        self.backscattering_coeff = np.zeros((0, 0, 0, 0))
        self.forwardscattering_coeff = np.zeros((0, 0, 0, 0))

        grids = [self.zenith_angles_scattering]

    def get_scattering_angles(self):
        return np.pi * self.zenith_angles_scattering / 180.0

    def get_absorption_coefficient(self):
        return self.extinction_matrix[0, 0, 0, 0, 0]

    def get_scattering_coefficient(self):
        return self.extinction_matrix[0, 0, 0, 0, 0] - self.absorption_vector[0, 0, 0, 0]

    def get_phase_function_integral(self):
        pf = self.phase_matrix[0, 0, 0, 0, :]
        x = np.cos(self.zenith_angles_scattering * np.pi / 180.0)

        norm = np.trapz(-2 * np.pi * pf, x=x)
        y = np.zeros(x.size)
        y[1:] = sp.integrate.cumtrapz(-2 * np.pi * pf, x=x)
        return y / norm

    def get_phase_function(self):
        pf = self.phase_matrix[0, 0, 0, 0, :]
        x = np.cos(self.zenith_angles_scattering * np.pi / 180.0)
        norm = np.trapz(-2 * np.pi * pf, x=x)
        return pf / norm

def setup_domain(n_cells):
    x = np.linspace(0, 100, n_cells + 1)
    y = np.linspace(-50, 50, n_cells + 1)
    z = np.linspace(-50, 50, n_cells + 1)
    return x, y, z

def get_absorption_coeff_field(n_particles, scattering_data):
    return n_particles * scattering_data.get_absorption_coefficient()

def get_scattering_coeff_field(n_particles, scattering_data):
    return n_particles * scattering_data.get_scattering_coefficient()

def get_phase_function_field(n_particles, scattering_data):
    n_dims = len(n_particles.shape)

    pf = scattering_data.get_phase_function_integral()
    return pf

scattering_data = ScatteringData("data/scattering_data_azimuthally_random.nc")
n = 201
x, y, z = setup_domain(n)
particles = np.zeros([n] * 3)
particles[n // 2, n // 2, n // 2] = 1e7
scattering_angles = scattering_data.get_scattering_angles()
absorption_coeff = get_absorption_coeff_field(particles, scattering_data)
scattering_coeff = get_scattering_coeff_field(particles, scattering_data)
phase_function = get_phase_function_field(particles, scattering_data)


ps = PointScattering(x, y, z, scattering_angles, absorption_coeff, scattering_coeff, phase_function)
in_pos, in_dir, out_pos, out_dir, scattering_events = ps.run(1000000)
print(np.sum(scattering_events==1))
out_dir = out_dir[scattering_events==1, :]
out_pos = out_pos[scattering_events==1, :]
theta= np.arccos(out_dir[:, 0])

bins = np.linspace(0, np.pi, 41)
y, _ = np.histogram(theta, bins=bins, density=True)
x = 0.5 * (bins[1:] + bins[:-1])
x_ref = scattering_data.get_scattering_angles()
y_ref = scattering_data.get_phase_function()
norm_ref = np.trapz(y_ref * np.sin(x_ref), x = x_ref)
plt.plot(x, y)
plt.plot(x_ref, y_ref * np.sin(x_ref) / norm_ref)
plt.show()

