import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

start_time = time.time()

dir1 = '/home/lex/Documents/Hexa/log/hexaRelax_2d_A_2048x2048_ex_ey_ne_0/run1'
dir2 = '/home/lex/Documents/Hexa/log/hexaRelax_2d_A_2048x2048_ex=ey=0/run1'
dir_list = [dir1, dir2]
n_dirs = len(dir_list)
dir_strings = [r'$E_x,\ E_y\ \ne\ 0$', '$E_x = E_y = 0$']
output_dir = 'figures/Video/c_along_x'

file_number = len(glob.glob(dir1 + '/relax_*'))

nx = 2048
nz = 2048
dx = np.float64(6 / nx)
dz = np.float64(6 / nz)
kx = np.pi / 3
Lz = 6
dt = 1e-2

iz_list = range(4)

var_list = ["ccx", "ccy", "ccz"]
for var in var_list:
    os.makedirs(output_dir + '/' + var, exist_ok = True)

cc = [0, 0]

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 4
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)

for n in range(file_number):

    dir_no = -1
    for dir in dir_list:
        dir_no += 1
        filename = dir + '/relax_' + '{:05d}'.format(n)
        file = FortranFile(filename, 'r')
        t = file.read_reals('float64')[0]
        bbx = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
        bby = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
        bbz = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")
        print(n, t)

        ccx = -(bby[1:, :] - bby[:-1, :]) / dz
        ccy = (bbx[1:, :] - bbx[:-1, :]) / dz - (bbz[:, 1:] - bbz[:, :-1]) / dx
        ccz = (bby[:, 1:] - bby[:, :-1]) / dx
        cc[dir_no] = {"ccx" : ccx, \
                      "ccy" : ccy, \
                      "ccz" : ccz}

    for var in var_list:
        for iz in iz_list:
            for dir_no in range(n_dirs):
                ax = fig.add_subplot(2, 4, iz + dir_no * 4 + 1)
                ax.plot(cc[dir_no][var][iz, :], \
                        label = dir_strings[dir_no])
                ax.set_title(var + ', iz = ' + str(iz))
                if iz == 1 and dir_no == 0:
                    ax.text(0.9, 1.1, \
                        	"t = " + "{:10.2e}".format(t), \
                            fontsize = 16, \
                        	transform = ax.transAxes)
                if iz == 3:
                    ax.text(1.1, 0.5, \
                        	dir_strings[dir_no], \
                            fontsize = 16, \
                        	transform = ax.transAxes)
        fig.savefig(output_dir + '/' + var + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
        fig.clf()

print(file_number)
print("--- %s seconds ---" % (time.time() - start_time))
