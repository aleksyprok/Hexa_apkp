import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

def ana_soln(x, z, var, l):
    if var == 'bbx':
        return -(l / kx) * np.cos(kx * x) * np.exp(-l * z)

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/along_z'

nx = 256
nz = 256
dx = np.float64(6 / nx)
dz = np.float64(6 / nz)
kx = np.pi / 3

xc = np.linspace(-dx / 2, 6 + dx / 2, nx + 2)
xb = np.linspace(0, 6, nx + 1)
zc = np.linspace(-dz / 2, 6 + dz / 2, nz + 2)
zb = np.linspace(-dz / 2, 6 + dz / 2, nz + 2)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ix_list = [0, 1, nx//4, nx//2, -nx//4, -2,-1]
var_list = ["bbx", "bby", "bbz", "divb"]
for var in var_list:
    os.makedirs(output_dir + '/' + var, exist_ok = True)

for n in range(0, file_number, 100):

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    t = file.read_reals('float64')[0]
    bbx = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
    bby = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
    bbz = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")

    print(n, t)

    divb = (bbx[1:-1, 1:] - bbx[1:-1, :-1]) / dx \
         + (bbz[1:, 1:-1] - bbz[:-1, 1:-1]) / dz

    bb = {"bbx" : bbx, \
          "bby" : bby, \
          "bbz" : bbz, \
          "divb" : divb}

    for var in var_list:
        k = 0
        for i in range(7):
            ix = ix_list[i]
            k += 1
            if i == 3: k += 1
            if i == 4: k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(bb[var][:, ix])
            ax.set_title(var + ', ix = ' + str(ix))
            if i == 1:
                ax.text(0.35, 1.1, \
                    	"t = " + "{:10.2e}".format(t), \
                        fontsize = 12, \
                    	transform=ax.transAxes)
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()
print("--- %s seconds ---" % (time.time() - start_time))
