import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os

def ay_ana_fun(x, z, l):
    return np.sin(kx * x) * np.sinh(l * (z - Lz)) / kx / np.sinh(-l * Lz)

def az_ana_fun(x, z, l):
    return np.sqrt(1 - l * l / (kx * kx)) * \
           np.cos(kx * x) * np.sinh(l * (z - Lz)) / kx / np.sinh(-l * Lz)

nx = 512
nz = 512
dx = np.float64(6 / nx)
dz = np.float64(6 / (nz + 1))
kx = np.pi / 3
Lz = 6

xc = np.linspace(-dx / 2, 6 + dx / 2, nx + 2)
xb = np.linspace(0, 6, nx + 1)
zc = np.linspace(dz / 2, 6 + dz / 2, nz + 2)
zb = np.linspace(dz, 6, nz + 1)

n = 42

filename = 'run1/relax_' + '{:05d}'.format(n)
file = FortranFile(filename, 'r')
t = file.read_reals('float64')[0]
bbx = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
bby = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
bbz = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")

bbx0 = bbx[0, :]
bbz1 = bbz[0, 1:-1]

alpha = 0.5 * kx
l = np.sqrt(kx * kx - alpha * alpha)
aay0 = ay_ana_fun(xb, 0, l)
bbz0 = (aay0[1:] - aay0[:-1]) / dx

divb0 = (bbx0[1:] - bbx0[:-1]) / dx \
      + (bbz1     - bbz0     ) / dz

fig = plt.figure()
ax = fig.add_subplot(111)
ln = ax.plot(divb0)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ln = ax.plot(bbz0)
ln = ax.plot(bbz1)
plt.show()
