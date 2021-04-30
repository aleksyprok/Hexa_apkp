import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

nx = 128
ny = 128
nz = 128
dx = np.float32(6 / nx)
dy = np.float32(6 / ny)
dz = np.float32(6 / nz)

filename = 'run1/relax_00000'
file = FortranFile(filename, 'r')
bbx0 = file.read_reals('float32').reshape((nz + 2, ny + 2, nx + 1), order = "C")
bby0 = file.read_reals('float32').reshape((nz + 2, ny + 1, nx + 2), order = "C")
bbz0 = file.read_reals('float32').reshape((nz + 1, ny + 2, nx + 2), order = "C")


filename = 'run1/relax_00999'
file = FortranFile(filename, 'r')
bbx = file.read_reals('float32').reshape((nz + 2, ny + 2, nx + 1), order = "C")
bby = file.read_reals('float32').reshape((nz + 2, ny + 1, nx + 2), order = "C")
bbz = file.read_reals('float32').reshape((nz + 1, ny + 2, nx + 2), order = "C")

ccx = (bbz[:, 1:, :] - bbz[:, :-1, :]) / dy \
    - (bby[1:, :, :] - bby[:-1, :, :]) / dz
ccy = (bbx[1:, :, :] - bbx[:-1, :, :]) / dz \
    - (bbz[:, :, 1:] - bbz[:, :, :-1]) / dx
ccz = (bby[:, :, 1:] - bby[:, :, :-1]) / dx \
    - (bbx[:, 1:, :] - bbx[:, :-1, :]) / dy

divb  = (bbx[1:-1, 1:-1, 1:] - bbx[1:-1, 1:-1, :-1]) / dx
divb += (bby[1:-1, 1:, 1:-1] - bby[1:-1, :-1, 1:-1]) / dy
divb += (bbz[1:, 1:-1, 1:-1] - bbz[:-1, 1:-1, 1:-1]) / dz

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(331)
ln = ax.plot(bbx[:, 32, 32], '-+')
ln = ax.plot(bbx0[:, 32, 32])

ax = fig.add_subplot(332)
ln = ax.plot(bbx[:, 32, 64], '-+')
ln = ax.plot(bbx0[:, 32, 64])

ax = fig.add_subplot(333)
ln = ax.plot(bbx[:, 32, 96], '-+')
ln = ax.plot(bbx0[:, 32, 96])

ax = fig.add_subplot(334)
ln = ax.plot(bbx[:, 64, 32], '-+')
ln = ax.plot(bbx0[:, 64, 32])

ax = fig.add_subplot(335)
ln = ax.plot(bbx[:, 64, 64], '-+')
ln = ax.plot(bbx0[:, 64, 64])

ax = fig.add_subplot(336)
ln = ax.plot(bbx[:, 64, 96], '-+')
ln = ax.plot(bbx0[:, 64, 96])

ax = fig.add_subplot(337)
ln = ax.plot(bbx[:, 96, 32], '-+')
ln = ax.plot(bbx0[:, 96, 32])

ax = fig.add_subplot(338)
ln = ax.plot(bbx[:, 96, 64], '-+')
ln = ax.plot(bbx0[:, 96, 64])

ax = fig.add_subplot(339)
ln = ax.plot(bbx[:, 96, 96], '-+')
ln = ax.plot(bbx0[:, 96, 96])

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(331)
ln = ax.plot(bby[:, 32, 32], '-+')
ln = ax.plot(bby0[:, 32, 32])

ax = fig.add_subplot(332)
ln = ax.plot(bby[:, 32, 64], '-+')
ln = ax.plot(bby0[:, 32, 64])

ax = fig.add_subplot(333)
ln = ax.plot(bby[:, 32, 96], '-+')
ln = ax.plot(bby0[:, 32, 96])

ax = fig.add_subplot(334)
ln = ax.plot(bby[:, 64, 32], '-+')
ln = ax.plot(bby0[:, 64, 32])

ax = fig.add_subplot(335)
ln = ax.plot(bby[:, 64, 64], '-+')
ln = ax.plot(bby0[:, 64, 64])

ax = fig.add_subplot(336)
ln = ax.plot(bby[:, 64, 96], '-+')
ln = ax.plot(bby0[:, 64, 96])

ax = fig.add_subplot(337)
ln = ax.plot(bby[:, 96, 32], '-+')
ln = ax.plot(bby0[:, 96, 32])

ax = fig.add_subplot(338)
ln = ax.plot(bby[:, 96, 64], '-+')
ln = ax.plot(bby0[:, 96, 64])

ax = fig.add_subplot(339)
ln = ax.plot(bby[:, 96, 96], '-+')
ln = ax.plot(bby0[:, 96, 96])

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(331)
ln = ax.plot(bbz[:, 32, 32], '-+')
ln = ax.plot(bbz0[:, 32, 32])

ax = fig.add_subplot(332)
ln = ax.plot(bbz[:, 32, 64], '-+')
ln = ax.plot(bbz0[:, 32, 64])

ax = fig.add_subplot(333)
ln = ax.plot(bbz[:, 32, 96], '-+')
ln = ax.plot(bbz0[:, 32, 96])

ax = fig.add_subplot(334)
ln = ax.plot(bbz[:, 64, 32], '-+')
ln = ax.plot(bbz0[:, 64, 32])

ax = fig.add_subplot(335)
ln = ax.plot(bbz[:, 64, 64], '-+')
ln = ax.plot(bbz0[:, 64, 64])

ax = fig.add_subplot(336)
ln = ax.plot(bbz[:, 64, 96], '-+')
ln = ax.plot(bbz0[:, 64, 96])

ax = fig.add_subplot(337)
ln = ax.plot(bbz[:, 96, 32], '-+')
ln = ax.plot(bbz0[:, 96, 32])

ax = fig.add_subplot(338)
ln = ax.plot(bbz[:, 96, 64], '-+')
ln = ax.plot(bbz0[:, 96, 64])

ax = fig.add_subplot(339)
ln = ax.plot(bbz[:, 96, 96], '-+')
ln = ax.plot(bbz0[:, 96, 96])

print(np.amax(np.abs(bbz[0, :, :])))
print(bbz[:, 64, 64])

plt.show()
