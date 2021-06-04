import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

def ramp_up(t):
    return np.sin(0.5 * np.pi * t) ** 2 * (t <= 1) + \
                                      1 * (t >  1)

def ax_ana_fun(x, z, l):
    return x - x

def ay_ana_fun(x, z, l):
    return np.sin(kx * x) * np.sinh(l * (z - Lz)) / kx / np.sinh(-l * Lz)

def az_ana_fun(x, z, l):
    return np.sqrt(1 - l * l / (kx * kx)) * np.cos(kx * x) * np.sinh(l * (z - Lz)) / kx / np.sinh(-l * Lz)

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/a_along_z'

nx = 256
nz = 256
dx = np.float64(6 / nx)
dz = np.float64(6 / nz)
kx = np.pi / 3
Lz = 6

xc = np.linspace(-dx / 2, 6 + dx / 2, nx + 2)
xb = np.linspace(0, 6, nx + 1)
zc = np.linspace(-dz / 2, 6 + dz / 2, nz + 2)
zb = np.linspace(0, 6, nz + 1)

alpha = 0.5 * kx
l = np.sqrt(kx * kx - alpha * alpha)

X, Z = np.meshgrid(xb, zc)
ax_ana = ax_ana_fun(X, Z, l)

X, Z = np.meshgrid(xc, zc)
ay_ana = ay_ana_fun(X, Z, l)

X, Z = np.meshgrid(xc, zb)
az_ana = az_ana_fun(X, Z, l)

a_ana = {"aax" : ax_ana, \
         "aay" : ay_ana, \
         "aaz" : az_ana}

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ix_list = [0, 1, nx//4, nx//2, -nx//4, -2,-1]
var_list = ["aax", "aay", "aaz"]
for var in var_list:
    os.makedirs(output_dir + '/' + var, exist_ok = True)

for n in range(0, file_number, 100):

    filename = 'run1/run1_' + '{:05d}'.format(n) + 'p'
    file = FortranFile(filename, 'r')
    opt = file.read_ints('int32')
    aax = file.read_reals('float32').reshape((nz + 1, nx), order = "C")
    aay = file.read_reals('float32').reshape((nz + 1, nx + 1), order = "C")
    aaz = file.read_reals('float32').reshape((nz, nx + 1), order = "C")

    print(n)

    aa = {"aax" : aax, \
          "aay" : aay, \
          "aaz" : aaz}

    for var in var_list:
        k = 0
        for i in range(7):
            ix = ix_list[i]
            k += 1
            if i == 3: k += 1
            if i == 4: k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(aa[var][:, ix])
            ax.plot(a_ana[var][:, ix])
            ax.set_title(var + ', ix = ' + str(ix))
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()
print("--- %s seconds ---" % (time.time() - start_time))
