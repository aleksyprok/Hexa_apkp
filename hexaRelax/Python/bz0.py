import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time

start_time = time.time()

nx = 128
ny = 128
nz = 128
dx = np.float32(6 / nx)
dy = np.float32(6 / ny)
dz = np.float32(6 / nz)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 1
fig.set_size_inches(fig_size)

filename = 'run1/poten_00002p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")
bbzp = (aay[:, :, 1:] - aay[:, :, :-1]) / dx \
     - (aax[:, 1:, :] - aax[:, :-1, :]) / dy
bbzp = bbzp[0, :, :]

filename = 'run1/run1_00001p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")
bbz1 = (aay[:, :, 1:] - aay[:, :, :-1]) / dx \
     - (aax[:, 1:, :] - aax[:, :-1, :]) / dy
bbz1 = bbz1[0, :, :]

filename = 'run1/run1_00002p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")
bbz2 = (aay[:, :, 1:] - aay[:, :, :-1]) / dx \
     - (aax[:, 1:, :] - aax[:, :-1, :]) / dy
bbz2 = bbz2[0, :, :]

filename = 'run1/run1_00003p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")
bbz3 = (aay[:, :, 1:] - aay[:, :, :-1]) / dx \
     - (aax[:, 1:, :] - aax[:, :-1, :]) / dy
bbz3 = bbz3[0, :, :]

ax = fig.add_subplot(131)
im = ax.imshow(bbz1 - bbzp, origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(132)
im = ax.imshow(bbz2 - bbzp, origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(133)
im = ax.imshow(bbz3 - bbzp, origin='lower')
cb = fig.colorbar(im)

print("--- %s seconds ---" % (time.time() - start_time))

plt.show()
