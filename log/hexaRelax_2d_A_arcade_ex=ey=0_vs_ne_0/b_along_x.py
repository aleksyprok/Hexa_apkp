import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

start_time = time.time()

dir1 = '/home/lex/Documents/Hexa/log/hexaRelax_2d_A_arcade_ex_ey_ne_0_first_100_iters/run1'
dir2 = '/home/lex/Documents/Hexa/log/hexaRelax_2d_A_arcade_ex=ey=0_first_100_iters/run1'
dir_list = [dir1, dir2]
n_dirs = len(dir_list)
dir_strings = ['ex_ey_ne_0', 'ex=ey=0']
output_dir = 'figures/Video/b_along_x'

file_number = len(glob.glob(dir1 + '/relax_*'))

nx = 256
nz = 256
dx = np.float64(6 / nx)
dz = np.float64(6 / nz)
kx = np.pi / 3
Lz = 6

iz_list = range(9)

var_list = ["bbx", "bby", "bbz", "b_mag"]
for var in var_list:
    os.makedirs(output_dir + '/' + var, exist_ok = True)

bb = [0, 0]

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

for n in range(3):

    dir_no = -1
    for dir in dir_list:
        dir_no += 1
        filename = dir + '/relax_' + '{:05d}'.format(n)
        file = FortranFile(filename, 'r')
        t = file.read_reals('float64')[0]
        bbx = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
        bby = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
        bbz = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")
        bx = 0.5 * (bbx[1:-1, :-1] + bbx[1:-1, 1:])
        by = bby[1:-1, 1:-1]
        bz = 0.5 * (bbz[:-1, 1:-1] + bbz[1:, 1:-1])
        b_mag = np.sqrt(bx * bx + by * by + bz * bz)
        print(n, t)
        bb[dir_no] = {"bbx" : bbx, \
                      "bby" : bby, \
                      "bbz" : bbz, \
                      "b_mag" : b_mag}

    for var in var_list:
        for iz in iz_list:
            ax = fig.add_subplot(3, 3, iz + 1)
            for dir_no in range(n_dirs):
                ax.plot(bb[1][var][iz, :], label = dir_strings[dir_no])
            ax.set_title(var + ', iz = ' + str(iz))
            if iz == 1:
                ax.text(0.35, 1.1, \
                    	"t = " + "{:10.2e}".format(t), \
                        fontsize = 12, \
                    	transform = ax.transAxes)
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()

print(file_number)
print("--- %s seconds ---" % (time.time() - start_time))
