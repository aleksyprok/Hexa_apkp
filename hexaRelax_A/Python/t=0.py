import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))

output_dir = 'Figures/Video/along_z'
os.makedirs(output_dir, exist_ok = True)
os.makedirs(output_dir + '/aax', exist_ok = True)
os.makedirs(output_dir + '/aay', exist_ok = True)
os.makedirs(output_dir + '/aaz', exist_ok = True)

nx = 128
ny = 128
nz = 128
dx = np.float64(6 / nx)
dy = np.float64(6 / ny)
dz = np.float64(6 / nz)

ix_list = [32, 64, 96]
iy_list = [32, 64, 96]

filename = 'setup_files/run1_00003p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)
k = 0
for i in range(3):
    ix = ix_list[i]
    for j in range(3):
        iy = iy_list[j]
        k += 1
        ax = fig.add_subplot(3, 3, k)
        ax.plot(aax[:, iy, ix])
        ax.set_title('aax, ix = ' + str(ix) + ', iy = ' + str(iy))

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)
k = 0
for i in range(3):
    ix = ix_list[i]
    for j in range(3):
        iy = iy_list[i]
        k += 1
        ax = fig.add_subplot(3, 3, k)
        ax.plot(aay[:, iy, ix])
        ax.set_title('aay, ix = ' + str(ix) + ', iy = ' + str(iy))

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)
k = 0
for i in range(3):
    ix = ix_list[i]
    for j in range(3):
        iy = iy_list[j]
        k += 1
        ax = fig.add_subplot(3, 3, k)
        ax.plot(aaz[:, iy, ix])
        ax.set_title('aaz, ix = ' + str(ix) + ', iy = ' + str(iy))

print("--- %s seconds ---" % (time.time() - start_time))

plt.show()
