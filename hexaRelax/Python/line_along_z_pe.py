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
dx = np.float64(6 / nx)
dy = np.float64(6 / ny)
dz = np.float64(6 / nz)

bbx_p = np.zeros((nz+2, ny+2, nx+1))
bby_p = np.zeros((nz+2, ny+1, nx+2))
bbz_p = np.zeros((nz+1, ny+2, nx+2))
bbx_e = np.zeros((nz+2, ny+2, nx+1))
bby_e = np.zeros((nz+2, ny+1, nx+2))
bbz_e = np.zeros((nz+1, ny+2, nx+2))

filename = 'setup_files/poten_00003p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")
bbx_p[1:-1, 1:-1, :] = (aaz[ :, 1:, :] - aaz[  :, :-1, :]) / dy \
                     - (aay[1:,  :, :] - aay[:-1,   :, :]) / dz
bby_p[1:-1, :, 1:-1] = (aax[1:, :,  :] - aax[:-1, :,   :]) / dz \
                     - (aaz[ :, :, 1:] - aaz[  :, :, :-1]) / dx
bbz_p[:, 1:-1, 1:-1] = (aay[:,  :, 1:] - aay[:,   :, :-1]) / dx \
                     - (aax[:, 1:,  :] - aax[:, :-1,   :]) / dy
bby_p[:, :, 0] = bby_p[:, :, 1]
bbz_p[:, :, 0] = bbz_p[:, :, 1]
bby_p[:, :, -1] = bby_p[:, :, -2]
bbz_p[:, :, -1] = bbz_p[:, :, -2]
bbx_p[:, 0, :] = bbx_p[:, 1, :]
bbz_p[:, 0, :] = bbz_p[:, 1, :]
bbx_p[:, -1, :] = bbx_p[:, -2, :]
bbz_p[:, -1, :] = bbz_p[:, -2, :]
bbx_p[0, :, :] = bbx_p[1, :, :] - dz * (bbz_p[0, :, 1:] - bbz_p[0, :, :-1]) / dx
bby_p[0, :, :] = bby_p[1, :, :] - dz * (bbz_p[0, 1:, :] - bbz_p[0, :-1, :]) / dy
bbx_p[-1, :, :] = bbx_p[-2, :, :]
bby_p[-1, :, :] = bby_p[-2, :, :]

filename = 'setup_files/run1_00003p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")
bbx_e[1:-1, 1:-1, :] = (aaz[ :, 1:, :] - aaz[  :, :-1, :]) / dy \
                     - (aay[1:,  :, :] - aay[:-1,   :, :]) / dz
bby_e[1:-1, :, 1:-1] = (aax[1:, :,  :] - aax[:-1, :,   :]) / dz \
                     - (aaz[ :, :, 1:] - aaz[  :, :, :-1]) / dx
bbz_e[:, 1:-1, 1:-1] = (aay[:,  :, 1:] - aay[:,   :, :-1]) / dx \
                     - (aax[:, 1:,  :] - aax[:, :-1,   :]) / dy
bby_e[:, :, 0] = bby_e[:, :, 1]
bbz_e[:, :, 0] = bbz_e[:, :, 1]
bby_e[:, :, -1] = bby_e[:, :, -2]
bbz_e[:, :, -1] = bbz_e[:, :, -2]
bbx_e[:, 0, :] = bbx_e[:, 1, :]
bbz_e[:, 0, :] = bbz_e[:, 1, :]
bbx_e[:, -1, :] = bbx_e[:, -2, :]
bbz_e[:, -1, :] = bbz_e[:, -2, :]
bbx_e[0, :, :] = bbx_e[1, :, :] - dz * (bbz_e[0, :, 1:] - bbz_e[0, :, :-1]) / dx
bby_e[0, :, :] = bby_e[1, :, :] - dz * (bbz_e[0, 1:, :] - bbz_e[0, :-1, :]) / dy
bbx_e[-1, :, :] = bbx_e[-2, :, :]
bby_e[-1, :, :] = bby_e[-2, :, :]

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
    bbx = file.read_reals('float64').reshape((nz + 2, ny + 2, nx + 1), order = "C")
    bby = file.read_reals('float64').reshape((nz + 2, ny + 1, nx + 2), order = "C")
    bbz = file.read_reals('float64').reshape((nz + 1, ny + 2, nx + 2), order = "C")

    k = 0
    for i in range(3):
        ix = ix_list[i]
        for j in range(3):
            iy = iy_list[j]
            k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(bbx[:, iy, ix])
            ax.plot(bbx_p[:, iy, ix])
            ax.plot(bbx_e[:, iy, ix])
            ax.set_xlim(0, 10)
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
            ax.plot(bby_p[:, iy, ix])
            ax.plot(bby_e[:, iy, ix])
            ax.set_xlim(0, 10)
            ax.set_title('bby, ix = ' + str(ix) + ', iy = ' + str(iy))
    fig.savefig(output_dir + '/bby/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()

    k = 0
    for i in range(3):
        ix = ix_list[i]
        for j in range(3):
            iy = iy_list[j]
            k += 1
            ax = fig.add_subplot(3, 3, k)
            ax.plot(bbz[:, iy, ix])
            ax.plot(bbz_p[:, iy, ix])
            ax.plot(bbz_e[:, iy, ix])
            ax.set_xlim(0, 10)
            ax.set_title('bbz, ix = ' + str(ix) + ', iy = ' + str(iy))
    fig.savefig(output_dir + '/bbz/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()

print("--- %s seconds ---" % (time.time() - start_time))
