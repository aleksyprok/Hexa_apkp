import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from streamtracer import StreamTracer, VectorGrid

nx = 128
x_min = -1
x_max = 1
dx = (x_max - x_min) / nx
xb = np.linspace(x_min, x_max, nx + 1)
xc = np.linspace(x_min - dx / 2, x_max + dx / 2, nx + 2)

ny = 128
y_min = -1
y_max = 1
dy = (y_max - y_min) / ny
yb = np.linspace(y_min, y_max, ny + 1)
yc = np.linspace(y_min - dy / 2, y_max + dy / 2, ny + 2)

nz = 128
z_min = 0
z_max = 2
dz = (z_max - z_min) / nz
zb = np.linspace(z_min, z_max, nz + 1)
zc = np.linspace(z_min - dz / 2, z_max + dz / 2, nz + 2)

# filename = 'setup_files/low_and_low_nlff.dat'
filename = 'run1/run1_' + '{:05d}'.format(1000)
file = FortranFile(filename, 'r')
bbx = file.read_reals('float64').reshape((nz + 2, ny + 2, nx + 1), order = "C").T
bby = file.read_reals('float64').reshape((nz + 2, ny + 1, nx + 2), order = "C").T
bbz = file.read_reals('float64').reshape((nz + 1, ny + 2, nx + 2), order = "C").T

bx = 0.25 * (bbx[:, :-1, :-1] + bbx[:, 1:, :-1] + \
             bbx[:, :-1,  1:] + bbx[:, 1:,  1:])
by = 0.25 * (bby[:-1, :, :-1] + bby[1:, :, :-1] + \
             bby[:-1, :,  1:] + bby[1:, :,  1:])
bz = 0.25 * (bbz[:-1, :-1, :] + bbz[1:, :-1, :] + \
             bbz[:-1,  1:, :] + bbz[1:,  1:, :])

nsteps = 10000
step_size = 0.01
tracer = StreamTracer(nsteps, step_size)

field = np.zeros((nx+1, ny+1, nz+1, 3))
field[:, :, :, 0] = bx
field[:, :, :, 1] = by
field[:, :, :, 2] = bz
grid_spacing = [dx, dy, dz]
grid = VectorGrid(field, grid_spacing, \
                  origin_coord = [x_min, y_min, z_min])

seeds = np.array([[-0.5,  0.0, 0.0], \
                  [ 0.0,  0.0, 0.0], \
                  [ 0.5,  0.0, 0.0], \
                  [ 0.0, -0.5, 0.0], \
                  [ 0.0,  0.5, 0.0], \
                  [-0.5,  0.0, 0.5], \
                  [ 0.0,  0.0, 0.5], \
                  [ 0.5,  0.0, 0.5], \
                  [ 0.0, -0.5, 0.5], \
                  [ 0.0,  0.5, 0.5]])
# seeds = np.array([[-0.5,  0.0, 1]])
tracer.trace(seeds, grid)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
X, Y = np.meshgrid(xc, yc)
cs = ax.contourf(X, Y, bbz[:, :, 0], \
                 zdir = 'z', \
                 offset = z_min, \
                 zorder = 0, \
                 levels = 30)
for n in range(10):
    ln = ax.plot(tracer.xs[n][:, 0], tracer.xs[n][:, 1], tracer.xs[n][:, 2], \
                 'black', \
                 zorder = 100)
ax.set_xlim((x_min, x_max))
ax.set_ylim((y_min, y_max))
ax.set_zlim((z_min, z_max))

plt.show()
