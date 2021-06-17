import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/check_b'

nx = 2048
nz = 2048
dx = np.float64(6 / nx)
dz = np.float64(6 / nz)
kx = np.pi / 3
Lz = 6
dt = 1e-2

xc = np.linspace(-dx / 2, 6 + dx / 2, nx + 2)
xb = np.linspace(0, 6, nx + 1)
zc = np.linspace(-dz / 2, 6 + dz / 2, nz + 2)
zb = np.linspace(0, 6, nz + 1)

filename = 'run1/relax_' + '{:05d}'.format(0)
file = FortranFile(filename, 'r')
t = file.read_reals('float64')[0]
bx_check = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
by_check = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
bz_check = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")

b_check = {"bbx" : bx_check, \
           "bby" : by_check, \
           "bbz" : bz_check}

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
            ax.plot(bb[var][iz, :])
            ax.plot(b_check[var][iz, :])
            ax.set_title(var + ', iz = ' + str(iz))
            if i == 1:
                ax.text(0.35, 1.1, \
                    	"t = " + "{:10.2e}".format(t), \
                        fontsize = 12, \
                    	transform=ax.transAxes)
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()

    frc = 1.9797798916682445e-3

    print(n, t)

    ccx = -(bby[1:, :] - bby[:-1, :]) / dz
    ccy = (bbx[1:, :] - bbx[:-1, :]) / dz - (bbz[:, 1:] - bbz[:, :-1]) / dx
    ccz = (bby[:, 1:] - bby[:, :-1]) / dx

    bx = 0.5 * (bbx[:-1, :] + bbx[1:, :])
    by = 0.25 * (bby[:-1, :-1] + bby[1:, :-1] + \
                 bby[:-1, 1: ] + bby[1:, 1: ])
    bz = 0.5 * (bbz[:, :-1] + bbz[:, 1:])
    bb = bx * bx + by * by + bz * bz

    bbmax = np.ones_like(bb) * 1e-8 * np.amax(bb)
    bbm = bb * (bb >= bbmax) + bbmax * (bb < bbmax)

    cx = 0.5 * (ccx[:, :-1] + ccx[:, 1:])
    cy = ccy
    cz = 0.5 * (ccz[:-1, :] + ccz[1:, :])

    vx = frc * (cy * bz - cz * by) / bbm
    vy = frc * (cz * bx - cx * bz) / bbm
    vz = frc * (cx * by - cy * bx) / bbm

    ex = vz * by - vy * bz
    ey = vx * bz - vz * bx
    ez = vy * bx - vx * by

    eex = 0.5 * (ex[:, :-1] + ex[:, 1:])
    eey = ey
    eez = 0.5 * (ez[:-1, :] + ez[1:, :])

    bx_check[1:-1, :] = bx_check[1:-1, :] + dt * (eey[1:, :] - eey[:-1, :]) / dz
    by_check[1:-1, 1:-1] = by_check[1:-1, 1:-1] + dt * ( \
                           (eez[:, 1:] - eez[:, :-1]) / dx \
                         - (eex[1:, :] - eex[:-1, :]) / dz)
    bz_check[:, 1:-1] = bz_check[:, 1:-1] - dt * (eey[:, 1:] - eey[:, :-1]) / dx

    b_check = {"bbx" : bx_check, \
               "bby" : by_check, \
               "bbz" : bz_check}

print("--- %s seconds ---" % (time.time() - start_time))
