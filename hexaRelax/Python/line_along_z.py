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
os.makedirs(output_dir + '/bbx', exist_ok = True)
os.makedirs(output_dir + '/bby', exist_ok = True)
os.makedirs(output_dir + '/bbz', exist_ok = True)

nx = 128
ny = 128
nz = 128
dx = np.float32(6 / nx)
dy = np.float32(6 / ny)
dz = np.float32(6 / nz)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)

ix_list = [32, 64, 96]
iy_list = [32, 64, 96]

for n in range(file_number):

    print(n)

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    bbx = file.read_reals('float32').reshape((nz + 2, ny + 2, nx + 1), order = "C")
    bby = file.read_reals('float32').reshape((nz + 2, ny + 1, nx + 2), order = "C")
    bbz = file.read_reals('float32').reshape((nz + 1, ny + 2, nx + 2), order = "C")

    k = 0
    for i in range(3):
        ix = ix_list[i]
        for j in range(3):
            iy = iy_list[j]
            k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(bbx[:, iy, ix])
            ax.set_title('bbx, ix = ' + str(ix) + ', iy = ' + str(iy))
    fig.savefig(output_dir + '/bbx/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()

    k = 0
    for i in range(3):
        ix = ix_list[i]
        for j in range(3):
            iy = iy_list[i]
            k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(bby[:, iy, ix])
            ax.set_title('bby, ix = ' + str(ix) + ', iy = ' + str(iy))
    fig.savefig(output_dir + '/bby/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()

    k = 0
    for i in range(3):
        ix = ix_list[i]
        for j in range(3):
            iy = iy_list[i]
            k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(bbz[:, iy, ix])
            ax.set_title('bbz, ix = ' + str(ix) + ', iy = ' + str(iy))
    fig.savefig(output_dir + '/bbz/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()

print("--- %s seconds ---" % (time.time() - start_time))
