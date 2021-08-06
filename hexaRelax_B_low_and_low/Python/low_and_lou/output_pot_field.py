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

filename = 'setup_files/low_and_low_nlff.dat'
file = FortranFile(filename, 'r')
bbx = file.read_reals('float64').reshape((nz + 2, ny + 2, nx + 1), order = "C").T
bby = file.read_reals('float64').reshape((nz + 2, ny + 1, nx + 2), order = "C").T
bbz = file.read_reals('float64').reshape((nz + 1, ny + 2, nx + 2), order = "C").T

maxIter = int(1e2)

# Set initial condition
# phi = np.zeros((nx+2, ny+2, nz+2))
filename = 'setup_files/low_and_low_pot_phi.dat'
file = FortranFile(filename, 'r')
phi = file.read_reals('float64').reshape((nz + 2, ny + 2, nx + 2), order = "C").T

# Apply boundary conditions
phi[ 0, :, :] = -bbx[ 0, :, :] * dx + phi[ 1, :, :]
phi[-1, :, :] =  bbx[-1, :, :] * dx + phi[-2, :, :]
phi[:,  0, :] = -bby[:,  0, :] * dy + phi[:,  1, :]
phi[:, -1, :] =  bby[:, -1, :] * dy + phi[:, -2, :]
phi[:, :,  0] = -bbz[:, :,  0] * dz + phi[:, :,  1]
phi[:, :, -1] =  bbz[:, :, -1] * dz + phi[:, :, -2]

for iteration in range(maxIter):
    if iteration % 100 == 0:
        bbx = (phi[1:, :, :] - phi[:-1, :, :]) / dx
        bby = (phi[:, 1:, :] - phi[:, :-1, :]) / dy
        bbz = (phi[:, :, 1:] - phi[:, :, :-1]) / dz
        divb = (bbx[1:, 1:-1, 1:-1] - bbx[:-1, 1:-1, 1:-1]) / dx \
             + (bby[1:-1, 1:, 1:-1] - bby[1:-1, :-1, 1:-1]) / dy \
             + (bbz[1:-1, 1:-1, 1:] - bbz[1:-1, 1:-1, :-1]) / dz
        print(iteration, np.max(np.abs(phi)), np.max(np.abs(divb)))
        f = FortranFile('setup_files/low_and_low_pot_phi.dat', 'w')
        f.write_record(phi.T)
        f.close()
    phi[1:-1, 1:-1, 1:-1] = (phi[2:, 1:-1, 1:-1] + phi[:-2, 1:-1, 1:-1] + \
                             phi[1:-1, 2:, 1:-1] + phi[1:-1, :-2, 1:-1] + \
                             phi[1:-1, 1:-1, 2:] + phi[1:-1, 1:-1, :-2]) / 6
    # Apply boundary conditions
    phi[ 0, :, :] = -bbx[ 0, :, :] * dx + phi[ 1, :, :]
    phi[-1, :, :] =  bbx[-1, :, :] * dx + phi[-2, :, :]
    phi[:,  0, :] = -bby[:,  0, :] * dy + phi[:,  1, :]
    phi[:, -1, :] =  bby[:, -1, :] * dy + phi[:, -2, :]
    phi[:, :,  0] = -bbz[:, :,  0] * dz + phi[:, :,  1]
    phi[:, :, -1] =  bbz[:, :, -1] * dz + phi[:, :, -2]

bbx = (phi[1:, :, :] - phi[:-1, :, :]) / dx
bby = (phi[:, 1:, :] - phi[:, :-1, :]) / dy
bbz = (phi[:, :, 1:] - phi[:, :, :-1]) / dz

f = FortranFile('setup_files/low_and_low_pot.dat', 'w')
f.write_record(bbx.T)
f.write_record(bby.T)
f.write_record(bbz.T)
f.close()

divb = (bbx[1:, 1:-1, 1:-1] - bbx[:-1, 1:-1, 1:-1]) / dx \
     + (bby[1:-1, 1:, 1:-1] - bby[1:-1, :-1, 1:-1]) / dy \
     + (bbz[1:-1, 1:-1, 1:] - bbz[1:-1, 1:-1, :-1]) / dz

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

seeds = np.array([[-0.5,  -0.0, 0.5], \
                  [-0.4,  -0.0, 0.5], \
                  [-0.3,  -0.0, 0.5], \
                  [-0.2,  -0.0, 0.5], \
                  [-0.1,  -0.0, 0.5], \
                  [ 0.0,  -0.0, 0.5], \
                  [ 0.1,  -0.0, 0.5], \
                  [ 0.2,  -0.0, 0.5], \
                  [ 0.3,  -0.0, 0.5], \
                  [ 0.4,  -0.0, 0.5], \
                  [ 0.5,  -0.0, 0.5]])
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

bx_p = bx
by_p = by
bz_p = bz
bx_pc = 0.5 * (bbx[:-1, 1:-1, 2:-1] + bbx[1:, 1:-1, 2:-1])
by_pc = 0.5 * (bby[1:-1, :-1, 2:-1] + bby[1:-1, 1:, 2:-1])
bz_pc = 0.5 * (bbz[1:-1, 1:-1, 1:-1] + bbz[1:-1, 1:-1, 2:])

##############################################################################

filename = 'setup_files/low_and_low_nlff.dat'
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

seeds = np.array([[-0.5,  -0.0, 0.5], \
                  [-0.4,  -0.0, 0.5], \
                  [-0.3,  -0.0, 0.5], \
                  [-0.2,  -0.0, 0.5], \
                  [-0.1,  -0.0, 0.5], \
                  [ 0.0,  -0.0, 0.5], \
                  [ 0.1,  -0.0, 0.5], \
                  [ 0.2,  -0.0, 0.5], \
                  [ 0.3,  -0.0, 0.5], \
                  [ 0.4,  -0.0, 0.5], \
                  [ 0.5,  -0.0, 0.5]])
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

bx_c = 0.5 * (bbx[:-1, 1:-1, 2:-1] + bbx[1:, 1:-1, 2:-1])
by_c = 0.5 * (bby[1:-1, :-1, 2:-1] + bby[1:-1, 1:, 2:-1])
bz_c = 0.5 * (bbz[1:-1, 1:-1, 1:-1] + bbz[1:-1, 1:-1, 2:])

print(np.max(np.sqrt((bx_p - bx)**2 + (by_p - by)**2 + (bz_p - bz)**2)))
print(np.sum(bx_pc**2 + by_pc**2 + bz_pc**2) * dx * dy * dz / 2)
print(np.sum(bx_c**2 + by_c**2 + bz_c**2) * dx * dy * dz / 2)
print(np.sum(bx_pc**2 + by_pc**2 + bz_pc**2) * dx * dy * dz / 2 - \
      np.sum(bx_c**2 + by_c**2 + bz_c**2) * dx * dy * dz / 2)

plt.show()
