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
var_list = ["bbx", "bby", "bbz", "divb"]
# var_list = ["diva"]
iz_list = [0, 1, 2, 3, nz//4, nz //2, 3*nz//4, -1, -2]
for var in var_list:
    os.makedirs(output_dir + '/' + var, exist_ok = True)

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

    aa = {"bbx" : bbx, \
          "bby" : bby, \
          "bbz" : bbz, \
          "divb" : divb}

    for var in var_list:
        for k in range(9):
            ax = fig.add_subplot(3, 3, k + 1)
            im = ax.imshow(aa[var][iz_list[k], :, :], origin='lower')
            cb = fig.colorbar(im)
            ax.set_title(var + ', iz = ' + str(iz_list[k]))
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()
print("--- %s seconds ---" % (time.time() - start_time))
