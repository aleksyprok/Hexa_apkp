import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

def jz_ana_fun(x, z, l):
    return np.sqrt(kx * kx - l * l) * np.cos(kx * x) * \
           np.sinh((l * (z - Lz))) / np.sinh(-l * Lz)

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/c_along_x'

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

alpha = 0.5 * kx
l = np.sqrt(kx * kx - alpha * alpha)
X, Z = np.meshgrid(xc, zb)
jz_ana = jz_ana_fun(X, Z, l)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

iz_list = [0, 1, 2, 3, nz//4, nz//2, -nz//4, -2,-1]
var_list = ["ccx", "ccy", "ccz"]
for var in var_list:
    os.makedirs(output_dir + '/' + var, exist_ok = True)

for n in range(0, file_number, 10):

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    t = file.read_reals('float64')[0]
    bbx = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
    bby = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
    bbz = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")

    print(n, t)

    ccx =-(bby[1:, :] - bby[:-1, :]) / dz
    ccy = (bbx[1:, :] - bbx[:-1, :]) / dz \
        - (bbz[:, 1:] - bbz[:, :-1]) / dx
    ccz = (bby[:, 1:] - bby[:, :-1]) / dx

    cc = {"ccx" : ccx, \
          "ccy" : ccy, \
          "ccz" : ccz}

    for var in var_list:
        k = 0
        for i in range(9):
            iz = iz_list[i]
            k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(cc[var][iz, :])
            if var == "ccz":
                ax.plot(jz_ana[iz, :])
            ax.set_title(var + ', iz = ' + str(iz))
            if i == 1:
                ax.text(0.35, 1.1, \
                    	"t = " + "{:10.2e}".format(t), \
                        fontsize = 12, \
                    	transform=ax.transAxes)
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()
print("--- %s seconds ---" % (time.time() - start_time))
