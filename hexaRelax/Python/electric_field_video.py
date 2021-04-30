import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time

start_time = time.time()

output_dir = 'Figures/Video/electric_x_y_plane'
os.makedirs(output_dir, exist_ok = True)

nx = 128
ny = 128
nz = 128
dx = np.float32(6 / nx)
dy = np.float32(6 / ny)
dz = np.float32(6 / nz)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)

for n in range(149):

    print(n)

    filename = 'run1/electric_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    eex = file.read_reals('float32').reshape((nz + 1, ny + 1, nx    ), order = "C")
    eey = file.read_reals('float32').reshape((nz + 1, ny,     nx + 1), order = "C")
    eez = file.read_reals('float32').reshape((nz,     ny + 1, nx + 1), order = "C")

    ax = fig.add_subplot(331)
    im = ax.imshow(eex[0, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('eex, z=0')

    ax = fig.add_subplot(332)
    im = ax.imshow(eey[0, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('eey, z=0')

    ax = fig.add_subplot(333)
    im = ax.imshow(eez[0, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('eez, z=0.5')

    ax = fig.add_subplot(334)
    im = ax.imshow(eex[1, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('eex, z=1')

    ax = fig.add_subplot(335)
    im = ax.imshow(eey[1, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('eey, z=1')

    ax = fig.add_subplot(336)
    im = ax.imshow(eez[1, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('eez, z=1.5')

    ax = fig.add_subplot(337)
    im = ax.imshow(eex[2, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('eex, z=2')

    ax = fig.add_subplot(338)
    im = ax.imshow(eey[2, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('eey, z=2')

    ax = fig.add_subplot(339)
    im = ax.imshow(eez[2, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('eez, z=2.5')

    fig.savefig(output_dir + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()

# output_dir = 'Figures/Video/db_dt_x_y_plane'
# os.makedirs(output_dir, exist_ok = True)
#
# for n in range(148):
#
#     print(n)
#
#     filename = 'run1/electric_' + '{:05d}'.format(n)
#     file = FortranFile(filename, 'r')
#     eex = file.read_reals('float32').reshape((nz + 1, ny + 1, nx    ), order = "C")
#     eey = file.read_reals('float32').reshape((nz + 1, ny,     nx + 1), order = "C")
#     eez = file.read_reals('float32').reshape((nz,     ny + 1, nx + 1), order = "C")
#
#     dbx_dt = (eey[1:, :, :] - eey[:-1, :, :]) / dz \
#            - (eez[:, 1:, :] - eez[:, :-1, :]) / dy
#     dby_dt = (eez[:, :, 1:] - eez[:, :, :-1]) / dx \
#            - (eex[1:, :, :] - eex[:-1, :, :]) / dz
#     dbz_dt = (eex[:, 1:, :] - eex[:, :-1, :]) / dy \
#            - (eey[:, :, 1:] - eey[:, :, :-1]) / dx
#
#     ax = fig.add_subplot(331)
#     im = ax.imshow(eex[0, :, :], origin='lower')
#     cb = fig.colorbar(im)
#     ax.set_title('dbx/dt, z=0.5')
#
#     ax = fig.add_subplot(332)
#     im = ax.imshow(eey[0, :, :], origin='lower')
#     cb = fig.colorbar(im)
#     ax.set_title('dby/dt, z=0.5')
#
#     ax = fig.add_subplot(333)
#     im = ax.imshow(eez[0, :, :], origin='lower')
#     cb = fig.colorbar(im)
#     ax.set_title('dbz/dt, z=0')
#
#     ax = fig.add_subplot(334)
#     im = ax.imshow(eex[1, :, :], origin='lower')
#     cb = fig.colorbar(im)
#     ax.set_title('dbx/dt, z=1.5')
#
#     ax = fig.add_subplot(335)
#     im = ax.imshow(eey[1, :, :], origin='lower')
#     cb = fig.colorbar(im)
#     ax.set_title('dby/dt, z=1.5')
#
#     ax = fig.add_subplot(336)
#     im = ax.imshow(eez[1, :, :], origin='lower')
#     cb = fig.colorbar(im)
#     ax.set_title('dbz/dt, z=1')
#
#     ax = fig.add_subplot(337)
#     im = ax.imshow(eex[2, :, :], origin='lower')
#     cb = fig.colorbar(im)
#     ax.set_title('dbx/dt, z=2.5')
#
#     ax = fig.add_subplot(338)
#     im = ax.imshow(eey[2, :, :], origin='lower')
#     cb = fig.colorbar(im)
#     ax.set_title('dby/dt, z=2.5')
#
#     ax = fig.add_subplot(339)
#     im = ax.imshow(eez[2, :, :], origin='lower')
#     cb = fig.colorbar(im)
#     ax.set_title('dbz/dt, z=2')
#
#     fig.savefig(output_dir + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
#     fig.clf()
