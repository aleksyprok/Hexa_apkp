import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

def bz_ana_fun(x, z, l):
    return np.cos(kx * x) * np.sinh((l * (z - Lz))) / np.sinh(-l * Lz)

def ay_ana_fun(x, z, l):
    return np.sin(kx * x) * np.sinh(l * (z - Lz)) / kx / np.sinh(-l * Lz)

nx = 512
nz = 512
dx = np.float64(6 / nx)
dz = np.float64(6 / (nz + 1))
kx = np.pi / 3
Lz = 6

xc = np.linspace(-dx / 2, 6 + dx / 2, nx + 2)
xb = np.linspace(0, 6, nx + 1)
zc = np.linspace(dz / 2, 6 + dz / 2, nz + 2)
zb = np.linspace(dz, 6, nz + 1)

alpha = 0.5 * kx
l = np.sqrt(kx * kx - alpha * alpha)
bz1 = bz_ana_fun(xc[1:-1], dz, l)

a = np.ones(nx) / (dx * dx)
b = - 2 * np.ones(nx) / (dx * dx)
c = np.ones(nx) / (dx * dx)
b0 = b[0]
b[0] += b0
b[-1] += a[0] * a[0] / b0

ab = np.zeros((3, nx))
ab[0,:] = a
ab[1,:] = b
ab[2,:] = a

u = np.zeros(nx)
v = np.zeros(nx)
u[0] = -b0
u[-1] = a[0]
v[0] = 1
v[-1] = -a[0] / b0

y = solve_banded((1,1), ab, -bz1)
q = solve_banded((1,1), ab, u)

vTy = v[0] * y[0] + v[-1] * y[-1]
vTq = v[0] * q[0] + v[-1] * q[-1]
P = y - vTy / (1 + vTq) * q

aay = np.zeros(nx)
aay[:-1] = -(P[1:] - P[:-1]) / dx
aay[-1] = -(P[0] - P[-1]) / dx

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 2
fig_size[1] = fig_size[1] * 1
fig.set_size_inches(fig_size)
ax = fig.add_subplot(121)
ln = ax.plot(aay)
ln = ax.plot(ay_ana_fun(xc[1:-1], dz, l))
ax = fig.add_subplot(122)
ln = ax.plot(aay - ay_ana_fun(xc[1:-1], dz, l))
plt.show()
