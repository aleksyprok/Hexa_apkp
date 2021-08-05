import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/along_xy'

nx = 128
ny = 128
nz = 128
dx = np.float64(6 / nx)
dy = np.float64(6 / ny)
dz = np.float64(6 / nz)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 4
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)

var_list = ["bbx", "bby", "bbz", "divb"]
for var in var_list:
    os.makedirs(output_dir + '/' + var, exist_ok = True)
ix = nx // 2
iy = ny // 2
iz_list = [0, 1, nz//2, -1]

for n in range(0, file_number, 10):

    print(n)

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    bbx = file.read_reals('float64').reshape((nz + 2, ny + 2, nx + 1), order = "C")
    bby = file.read_reals('float64').reshape((nz + 2, ny + 1, nx + 2), order = "C")
    bbz = file.read_reals('float64').reshape((nz + 1, ny + 2, nx + 2), order = "C")

    divb = (bbx[1:-1, 1:-1, 1:] - bbx[1:-1, 1:-1, :-1]) / dx \
         + (bby[1:-1, 1:, 1:-1] - bby[1:-1, :-1, 1:-1]) / dy \
         + (bbz[1:, 1:-1, 1:-1] - bbz[:-1, 1:-1, 1:-1]) / dz

    bb = {"bbx" : bbx, \
          "bby" : bby, \
          "bbz" : bbz, \
          "divb" : divb}

    for var in var_list:
        k = 0
        for iz in iz_list:
            k += 1
            ax = fig.add_subplot(2, 4, k)
            ax.plot(bb[var][iz, iy, :])
            ax.set_title(var + ', iy = ' + str(iy) + ', iz = ' + str(iz))
        for iz in iz_list:
            k += 1
            ax = fig.add_subplot(2, 4, k)
            ax.plot(bb[var][iz, :, ix])
            ax.set_title(var + ', ix = ' + str(ix) + ', iz = ' + str(iz))
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()
print("--- %s seconds ---" % (time.time() - start_time))
