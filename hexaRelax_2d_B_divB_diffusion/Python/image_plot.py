import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/image_plot'
os.makedirs(output_dir, exist_ok = True)

nx = 256
nz = 256
dx = np.float64(6 / nx)
dz = np.float64(6 / nz)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 2
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)

var_list = ["bbx", "bby", "bbz", "divb"]

for n in range(0, file_number, 100):

    print(n)

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    t = file.read_reals('float64')[0]
    bbx = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
    bby = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
    bbz = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")

    divb = (bbx[1:-1, 1:] - bbx[1:-1, :-1]) / dx \
         + (bbz[1:, 1:-1] - bbz[:-1, 1:-1]) / dz

    bb = {"bbx" : bbx, \
          "bby" : bby, \
          "bbz" : bbz, \
          "divb" : divb}

    k = 0
    for var in var_list:
        k += 1
        ax = fig.add_subplot(2, 2, k)
        im = ax.imshow(bb[var], origin='lower')
        cb = fig.colorbar(im)
        ax.set_title(var)
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        if k == 1:
            ax.text(1.2, 1.1, \
                	"t = " + "{:10.2e}".format(t), \
                    fontsize = 12, \
                	transform=ax.transAxes)
    fig.savefig(output_dir + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()
print("--- %s seconds ---" % (time.time() - start_time))
