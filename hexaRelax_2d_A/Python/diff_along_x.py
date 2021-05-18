import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/diff_along_x'

nx = 2048
nz = 2048
dx = np.float64(6 / nx)
dz = np.float64(6 / nz)
kx = np.pi / 3
Lz = 6

xc = np.linspace(-dx / 2, 6 + dx / 2, nx + 2)
xb = np.linspace(0, 6, nx + 1)
zc = np.linspace(-dz / 2, 6 + dz / 2, nz + 2)
zb = np.linspace(0, 6, nz + 1)

filename = 'run1/relax_' + '{:05d}'.format(0)
file = FortranFile(filename, 'r')
t = file.read_reals('float64')[0]
bbx = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
bby = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
bbz = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")
bb0 = {"bbx" : bbx, \
       "bby" : bby, \
       "bbz" : bbz}

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

iz_list = [0, 1, nz//4, nz//2, -nz//4, -2, -1]
var_list = ["bbx", "bby", "bbz"]
for var in var_list:
    os.makedirs(output_dir + '/' + var, exist_ok = True)

for n in range(0, file_number, 1):

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    t = file.read_reals('float64')[0]
    bbx = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
    bby = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
    bbz = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")

    print(n, t)

    bb = {"bbx" : bbx, \
          "bby" : bby, \
          "bbz" : bbz}

    for var in var_list:
        k = 0
        for i in range(7):
            iz = iz_list[i]
            k += 1
            if i == 3: k += 1
            if i == 4: k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(bb[var][iz, :] - bb0[var][iz, :])
            ax.set_title(var + ' - ' + var + '0' + ', iz = ' + str(iz))
            if i == 1:
                ax.text(0.35, 1.1, \
                    	"t = " + "{:10.2e}".format(t), \
                        fontsize = 12, \
                    	transform=ax.transAxes)
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()
print("--- %s seconds ---" % (time.time() - start_time))
