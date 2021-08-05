import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp

def dr_ds(s, r):
    # r[0] = x
    # r[1] = y
    # r[2] = z
    # dx_ds = bx_interp([r[0], r[1], r[2]])
    # dy_ds = by_interp([r[0], r[1], r[2]])
    # dz_ds = bz_interp([r[0], r[1], r[2]])
    dx_ds = r[0]
    dy_ds = r[1]
    dz_ds = r[2]
    return np.array([dx_ds, dy_ds, dz_ds])

nx = 128
x_min = -0.8
x_max = 1.2
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
z_min = 1
z_max = 3
dz = (z_max - z_min) / nz
zb = np.linspace(z_min, z_max, nz + 1)
zc = np.linspace(z_min - dz / 2, z_max + dz / 2, nz + 2)

filename = 'run1/run1_' + '{:05d}'.format(0)
file = FortranFile(filename, 'r')
bbx = file.read_reals('float64').reshape((nz + 2, ny + 2, nx + 1), order = "C").T
bby = file.read_reals('float64').reshape((nz + 2, ny + 1, nx + 2), order = "C").T
bbz = file.read_reals('float64').reshape((nz + 1, ny + 2, nx + 2), order = "C").T

bx_interp = RegularGridInterpolator((xb, yc, zc), bbx)
by_interp = RegularGridInterpolator((xc, yb, zc), bby)
bz_interp = RegularGridInterpolator((xc, yc, zb), bbz)

s_range = [0, 10]
r0 = [0, 0, 0]
sol = solve_ivp(dr_ds, s_range, r0)

# print(bz_interp([xc[5], yc[9], zb[6]]))
# print(bbz[5, 9, 6])

# X, Y = np.meshgrid(xc, yc)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# im = ax.imshow(bz_interp(X, Y, z_min), \
#                extent = [x_min, x_max, y_min, y_max], \
#                origin = 'lower')
# cb = fig.colorbar(im)

# X, Y = np.meshgrid(xc, yc)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection = '3d')
# cs = ax.contourf(X, Y, bbz[:, :, 0], \
#                  zdir = 'z', \
#                  offset = z_min, \
#                  levels = 30)
# ax.set_zlim((z_min, z_max))

plt.show()
