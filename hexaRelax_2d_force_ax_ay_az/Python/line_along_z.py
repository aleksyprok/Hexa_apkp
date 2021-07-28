import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

def ramp_up(t):
    return np.sin(0.5 * np.pi * t) ** 2 * (t <= 1) + \
                                      1 * (t >  1)

def bx_ana_fun(x, z, l):
    return -(l / kx) * np.sin(kx * x) * np.cosh((l * (z - Lz))) / np.sinh(-l * Lz)

def by_ana_fun(x, z, l):
    return np.sqrt(1 - l * l / (kx * kx)) * np.sin(kx * x) * np.sinh((l * (z - Lz))) / np.sinh(-l * Lz)

def bz_ana_fun(x, z, l):
    return np.cos(kx * x) * np.sinh((l * (z - Lz))) / np.sinh(-l * Lz)

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/along_z'

nx = 512
nz = 512
dx = np.float64(6 / nx)
dz = np.float64(6 / (nz + 1))
kx = np.pi / 3
Lz = 6

xc = np.linspace(-dx / 2, 6 + dx / 2, nx + 2)
xb = np.linspace(0, 6, nx + 1)
zc = np.linspace(dz / 2, 6 + dz / 2, nz + 2)
zb = np.linspace(dz, 6, nz + 1)

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

    divb = (bbx[1:-1,1:] - bbx[1:-1,:-1]) / dx \
         + (bbz[1:,1:-1] - bbz[:-1,1:-1]) / dz

    print(n, t)

    bb = {"bbx" : bbx, \
          "bby" : bby, \
          "bbz" : bbz, \
          "divb" : divb}

    alpha = 0.5 * kx * ramp_up(0.1 * t)
    l = np.sqrt(kx * kx - alpha * alpha)

    X, Z = np.meshgrid(xb, zc)
    bx_ana = bx_ana_fun(X, Z, l)

    X, Z = np.meshgrid(xc, zc)
    by_ana = by_ana_fun(X, Z, l)

    X, Z = np.meshgrid(xc, zb)
    bz_ana = bz_ana_fun(X, Z, l)

    b_ana = {"bbx" : bx_ana, \
             "bby" : by_ana, \
             "bbz" : bz_ana}

    for var in var_list:
        k = 0
        for i in range(7):
            ix = ix_list[i]
            k += 1
            if i == 3: k += 1
            if i == 4: k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(bb[var][:, ix])
            # ax.plot(bb[var][:5, ix])
            if var != "divb":
                ax.plot(b_ana[var][:, ix])
                # ax.plot(b_ana[var][:5, ix])
            ax.set_title(var + ', ix = ' + str(ix))
            if i == 1:
                ax.text(0.35, 1.1, \
                    	"t = " + "{:10.2e}".format(t), \
                        fontsize = 12, \
                    	transform=ax.transAxes)
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()
print("--- %s seconds ---" % (time.time() - start_time))
