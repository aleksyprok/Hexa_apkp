import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os

os.makedirs('Figures', exist_ok = True)

nx = 128
ny = 128
nz = 128
dx = np.float32(6 / nx)
dy = np.float32(6 / ny)
dz = np.float32(6 / nz)

filename = 'run1/relax_00000'
file = FortranFile(filename, 'r')
bbx = file.read_reals('float32').reshape((nz + 2, ny + 2, nx + 1), order = "C")
bby = file.read_reals('float32').reshape((nz + 2, ny + 1, nx + 2), order = "C")
bbz = file.read_reals('float32').reshape((nz + 1, ny + 2, nx + 2), order = "C")

filename = 'run1/poten_00003p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")

bbx_a = (aaz[ :, 1:, :] - aaz[  :, :-1, :]) / dy \
      - (aay[1:,  :, :] - aay[:-1,   :, :]) / dz
bby_a = (aax[1:, :,  :] - aax[:-1, :,   :]) / dz \
      - (aaz[ :, :, 1:] - aaz[  :, :, :-1]) / dx
bbz_a = (aay[:,  :, 1:] - aay[0,   :, :-1]) / dx \
      - (aax[:, 1:,  :] - aax[0, :-1,   :]) / dy

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
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)

ax = fig.add_subplot(341)
im = ax.imshow(ccx[0, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('ccx, z=0')

ax = fig.add_subplot(342)
im = ax.imshow(ccy[0, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('ccy, z=0')

ax = fig.add_subplot(343)
im = ax.imshow(ccz[0, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('ccz, z=-1/2')

ax = fig.add_subplot(344)
im = ax.imshow(divb[0, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('div B, z=1/2')

ax = fig.add_subplot(345)
im = ax.imshow(ccx[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(346)
im = ax.imshow(ccy[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(347)
im = ax.imshow(ccz[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(348)
im = ax.imshow(divb[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(349)
im = ax.imshow(ccx[2, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(3, 4, 10)
im = ax.imshow(ccy[2, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(3, 4, 11)
im = ax.imshow(ccz[2, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(3, 4, 12)
im = ax.imshow(divb[2, :, :], origin='lower')
cb = fig.colorbar(im)

fig.savefig('Figures/read_data1.png',  bbox_inches='tight')

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(331)
im = ax.imshow(bbx[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(332)
im = ax.imshow(bby[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(333)
im = ax.imshow(bbz[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(334)
im = ax.imshow(bbx[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(335)
im = ax.imshow(bby[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(336)
im = ax.imshow(bbz[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(337)
im = ax.imshow(bbx[2, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(338)
im = ax.imshow(bby[2, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(339)
im = ax.imshow(bbz[2, :, :], origin='lower')
cb = fig.colorbar(im)

fig.savefig('Figures/read_data2.png',  bbox_inches='tight')

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(331)
im = ax.imshow(bbx[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(332)
im = ax.imshow(bby[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(333)
im = ax.imshow(bbz[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(334)
im = ax.imshow(bbx_a[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(335)
im = ax.imshow(bby_a[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(336)
im = ax.imshow(bbz_a[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(337)
im = ax.imshow(bbx[1, 1:-1, :] - bbx_a[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(338)
im = ax.imshow(bby[1, :, 1:-1] - bby_a[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(339)
im = ax.imshow(bbz[0, 1:-1, 1:-1] - bbz_a[0, :, :], origin='lower')
cb = fig.colorbar(im)

fig.savefig('Figures/read_data3.png',  bbox_inches='tight')

# print(bbx[1, 1:-1, :] - bbx_a[0, :, :])
# print(bby[1, :, 1:-1] - bby_a[0, :, :])
# print(bbz[0, 1:-1, 1:-1] - bbz_a[0, :, :])
#
# print(np.amax(np.abs(bbx[1, 1:-1, :] - bbx_a[0, :, :])))
# print(np.amax(np.abs(bby[1, :, 1:-1] - bby_a[0, :, :])))
# print(np.amax(np.abs(bbz[0, 1:-1, 1:-1] - bbz_a[0, :, :])))
