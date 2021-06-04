import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time
import glob

start_time = time.time()

file_number = len(glob.glob('run1/relax_*'))
print(file_number)

output_dir = 'Figures/Video/fieldlines'
os.makedirs(output_dir, exist_ok = True)

nx = 2048
nz = 2048
dx = np.float64(6 / nx)
dz = np.float64(6 / nz)

x_min = 0
x_max = 6
x = np.linspace(dx / 2, x_max - dx / 2, nx)

z_min = 0
z_max = 6
z = np.linspace(dz / 2, z_max - dz / 2, nz)

X, Z = np.meshgrid(x, z)

start_points = []
x_coords = np.array([1.5, 4.5])
nzs = 10
dzs = np.float64(z_max / nzs)
z_coords = np.linspace(dzs / 2, z_max - dzs/2, nzs)
for xs in x_coords:
    for zs in z_coords:
        start_points.append([xs, zs])
start_points = np.array(start_points)

fig = plt.figure()

for n in range(0, file_number, 1):

    ax = fig.add_subplot(111)
    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    t = file.read_reals('float64')[0]
    bbx = file.read_reals('float64').reshape((nz + 2, nx + 1), order = "C")
    bby = file.read_reals('float64').reshape((nz + 2, nx + 2), order = "C")
    bbz = file.read_reals('float64').reshape((nz + 1, nx + 2), order = "C")

    bx = 0.5 * (bbx[1:-1, :-1] + bbx[1:-1, 1:])
    bz = 0.5 * (bbz[:-1, 1:-1] + bbz[1:, 1:-1])

    strm = ax.streamplot(X, Z, bx, bz, start_points = start_points, density = 35)
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title('Field lines, t =' + "{:10.2e}".format(t))

    fig.savefig(output_dir + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()

print("--- %s seconds ---" % (time.time() - start_time))
