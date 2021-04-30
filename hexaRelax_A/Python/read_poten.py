import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob


start_time = time.time()
# read_dir = '/home/lex/OneDrive/Documents/Hexa/hexa_freeze_bz/potenMPI/run1/'
read_dir = '/home/lex/OneDrive/Documents/Hexa/hexa_freeze_bz/hexaf90/run1/'

nx = 128
ny = 128
nz = 128
dx = np.float64(6 / nx)
dy = np.float64(6 / ny)
dz = np.float64(6 / nz)

filename = read_dir + 'run1_00060p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")
diva = (aax[1:-1, 1:-1, 1:] - aax[1:-1, 1:-1, :-1]) / dx \
     + (aay[1:-1, 1:, 1:-1] - aay[1:-1, :-1, 1:-1]) / dy \
     + (aaz[1:, 1:-1, 1:-1] - aaz[:-1, 1:-1, 1:-1]) / dz

aa = {"aax" : aax, \
      "aay" : aay, \
      "aaz" : aaz, \
      "diva" : diva}

ix_list = [32, 64, 96]
iy_list = [32, 64, 96]
var_list = ["aax", "aay", "aaz", "diva"]

for var in var_list:
    fig = plt.figure()
    fig_size = fig.get_size_inches()
    fig_size[0] = fig_size[0] * 3
    fig_size[1] = fig_size[1] * 3
    fig.set_size_inches(fig_size)
    k = 0
    for i in range(3):
        ix = ix_list[i]
        for j in range(3):
            iy = iy_list[j]
            k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(aa[var][:, iy, ix])
            ax.set_title(var + ', ix = ' + str(ix) + ', iy = ' + str(iy))

plt.show()
