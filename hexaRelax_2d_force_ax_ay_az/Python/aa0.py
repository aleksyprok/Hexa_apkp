import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

file_number = len(glob.glob('run1/aa0_*'))
print(file_number)

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

n = 1

filename = 'run1/aa0_' + '{:05d}'.format(n)
file = FortranFile(filename, 'r')
t = file.read_reals('float64')[0]
aax0 = file.read_reals('float64')
aay0 = file.read_reals('float64')
aaz0 = file.read_reals('float64')
aax1 = file.read_reals('float64')
aay1 = file.read_reals('float64')

filename = 'run1/run1_' + '{:05d}'.format(n) + 'p'
file = FortranFile(filename, 'r')
opt = file.read_ints('int32')
aax = file.read_reals('float32').reshape((nz + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, nx + 1), order = "C")

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(331)
ax.plot(aax1)
ax = fig.add_subplot(332)
ax.plot(aay1)
ax = fig.add_subplot(333)
ax.plot(aaz0)

ax = fig.add_subplot(334)
ax.plot(aax1)
ax = fig.add_subplot(335)
ax.plot(aay1)

ax = fig.add_subplot(337)
ax.plot(aax[0, :])
ax = fig.add_subplot(338)
ax.plot(aay[0, :])

plt.show()
