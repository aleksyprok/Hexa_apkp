import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

nx = 256
ny = 256
nz = 256
dx = np.float32(6 / nx)
dy = np.float32(6 / ny)
dz = np.float32(6 / nz)

filename = 'setup_files/poten_00010p'
# filename = 'setup_files/128x128x128/poten_00003p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")

bbx_p = (aaz[ :, 1:, :] - aaz[  :, :-1, :]) / dy \
      - (aay[1:,  :, :] - aay[:-1,   :, :]) / dz
bby_p = (aax[1:, :,  :] - aax[:-1, :,   :]) / dz \
      - (aaz[ :, :, 1:] - aaz[  :, :, :-1]) / dx
bbz_p = (aay[:,  :, 1:] - aay[:,   :, :-1]) / dx \
      - (aax[:, 1:,  :] - aax[:, :-1,   :]) / dy

ccx_p = (bbz_p[1:-1, 1:, :] - bbz_p[1:-1, :-1, :]) / dy \
      - (bby_p[1:, 1:-1, :] - bby_p[:-1, 1:-1, :]) / dz
ccy_p = (bbx_p[1:, :, 1:-1] - bbx_p[:-1, :, 1:-1]) / dz \
      - (bbz_p[1:-1, :, 1:] - bbz_p[1:-1, :, :-1]) / dx
ccz_p = (bby_p[:, 1:-1, 1:] - bby_p[:, 1:-1, :-1]) / dx \
      - (bbx_p[:, 1:, 1:-1] - bbx_p[:, :-1, 1:-1]) / dy

filename = 'run1/relax_00000'
# filename = '/home/lex/OneDrive/Documents/Hexa/hexaRelax_potential_field_unstable/run1/relax_00000'
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


fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 2
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(321)
im = ax.imshow(bbx[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('bbx')

ax = fig.add_subplot(322)
im = ax.imshow(bbx_p[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('bbx_p')

ax = fig.add_subplot(323)
im = ax.imshow(bby[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('bby')

ax = fig.add_subplot(324)
im = ax.imshow(bby_p[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('bby_p')

ax = fig.add_subplot(325)
im = ax.imshow(bbz[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('bbz')

ax = fig.add_subplot(326)
im = ax.imshow(bbz_p[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('bbz_p')

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 2
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(321)
im = ax.imshow(ccx[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('ccx')

ax = fig.add_subplot(322)
im = ax.imshow(ccx_p[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('ccx_p')

ax = fig.add_subplot(323)
im = ax.imshow(ccy[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('ccy')

ax = fig.add_subplot(324)
im = ax.imshow(ccy_p[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('ccy_p')

ax = fig.add_subplot(325)
im = ax.imshow(ccz[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('ccz')

ax = fig.add_subplot(326)
im = ax.imshow(ccz_p[100, :, :], origin='lower')
cb = fig.colorbar(im)
ax.set_title('ccz_p')

plt.show()
