import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/contour_xy'

nx = 128
ny = 128
nz = 128
dx = np.float64(6 / nx)
dy = np.float64(6 / ny)
dz = np.float64(6 / nz)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ix_list = [32, 64, 96]
iy_list = [32, 64, 96]
var_list = ["bbx", "bby", "bbz"]
for var in var_list:
    os.makedirs(output_dir + '/' + var, exist_ok = True)

for n in range(0, file_number):

    print(n)

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    aax = file.read_reals('float64').reshape((nz + 1, ny + 1, nx), order = "C")
    aay = file.read_reals('float64').reshape((nz + 1, ny, nx + 1), order = "C")
    aaz = file.read_reals('float64').reshape((nz, ny + 1, nx + 1), order = "C")

    bbx = (aaz[:, 1:, :] - aaz[:, :-1, :]) / dy \
        - (aay[1:, :, :] - aay[:-1, :, :]) / dz
    bby = (aax[1:, :, :] - aax[:-1, :, :]) / dz \
        - (aaz[:, :, 1:] - aaz[:, :, :-1]) / dx
    bbz = (aay[:, :, 1:] - aay[:, :, :-1]) / dx \
        - (aax[:, 1:, :] - aax[:, :-1, :]) / dy

    bb = {"bbx" : bbx, \
          "bby" : bby, \
          "bbz" : bbz}

    for var in var_list:
        for k in range(9):
            ax = fig.add_subplot(3, 3, k + 1)
            im = ax.imshow(bb[var][k, :, :], origin='lower')
            cb = fig.colorbar(im)
            ax.set_title(var + ', iz = ' + str(k))
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()
print("--- %s seconds ---" % (time.time() - start_time))
