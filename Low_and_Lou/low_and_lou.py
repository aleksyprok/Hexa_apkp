import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

def dY_dmu(mu, Y, p):
    # Let Q = P'(mu)
    # Y[0] = P(mu)
    # Y[1] = Q(mu)
    n = 1
    a = p[0]
    dP_dmu = Y[1]
    dQ_dmu = n * (n + 1) * Y[0] + a * a * (1 + n) / n * Y[0] ** (1 + 2 / n)
    dQ_dmu = -dQ_dmu / (1 - mu * mu)
    return np.array([dP_dmu, dQ_dmu])

def bcs(Y_a, Y_b, p):
	return np.array([Y_a[0], Y_b[0], Y_a[1] - 10])

def B_r(r, theta):
    mu = np.cos(theta)
    return -sol.sol(mu)[1] / r ** 2

def B_theta(r, theta):
    mu = np.cos(theta)
    return sol.sol(mu)[0] / r ** 2

def B_phi(r, theta):
    mu = np.cos(theta)
    return -sol.p[0] * sol.sol(mu)[0] ** 2 / r ** 2

def B_X(X, Y, Z):
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    theta = np.arccos(Z / r)
    phi = np.arctan2(Y, X)
    return B_r(r, theta) * np.sin(theta) * np.cos(phi) + \
           B_theta(r, theta) * np.cos(theta) * np.cos(phi) - \
           B_phi(r, theta) * np.sin(phi)

def B_Y(X, Y, Z):
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    theta = np.arccos(Z / r)
    phi = np.arctan2(Y, X)
    return B_r(r, theta) * np.sin(theta) * np.sin(phi) + \
           B_theta(r, theta) * np.cos(theta) * np.sin(phi) + \
           B_phi(r, theta) * np.cos(phi)

def B_Z(X, Y, Z):
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    theta = np.arccos(Z / r)
    phi = np.arctan2(Y, X)
    return B_r(r, theta) * np.cos(theta) - \
           B_theta(r, theta) * np.sin(theta)

def B_x(x, y, z):
    X = (x + lx) * np.cos(Phi) - (z + lz) * np.sin(Phi)
    Y = y
    Z = (x + lx) * np.sin(Phi) + (z + lz) * np.cos(Phi)
    return B_X(X, Y, Z) * np.cos(Phi) + \
           B_Z(X, Y, Z) * np.sin(Phi)

def B_y(x, y, z):
    X = (x + lx) * np.cos(Phi) - (z + lz) * np.sin(Phi)
    Y = y
    Z = (x + lx) * np.sin(Phi) + (z + lz) * np.cos(Phi)
    return B_Y(X, Y, Z)

def B_z(x, y, z):
    X = (x + lx) * np.cos(Phi) - (z + lz) * np.sin(Phi)
    Y = y
    Z = (x + lx) * np.sin(Phi) + (z + lz) * np.cos(Phi)
    return -B_X(X, Y, Z) * np.sin(Phi) + \
            B_Z(X, Y, Z) * np.cos(Phi)


# Calculate P(mu) and P'(mu)
mu_min = -1 + 1e-12
mu_max = 1 - 1e-12
mu_temp = np.linspace(mu_min, mu_max, 5)
Y = np.zeros((2, mu_temp.size))
Y[0, 1] = 1
Y[0, 3] = -1
sol = solve_bvp(dY_dmu, bcs, mu_temp, Y, p = [np.sqrt(0.425)], tol = 1e-10)

lx = -0.25
lz = 0.0
Phi = np.pi / 10

nx = 128
x_min = 1
x_max = 3
x = np.linspace(x_min, x_max, nx)

ny = 128
y_min = -1
y_max = 1
y = np.linspace(y_min, y_max, ny)

nz = 128
z_min = -0.8
z_max = 1.2
z = np.linspace(z_min, z_max, nz)

n_mu = 1024
mu = np.linspace(mu_min, mu_max, n_mu)

x_grid, y_grid = np.meshgrid(x, y)
fig = plt.figure()
ax = fig.add_subplot(111)
# ct = ax.contour(x_grid, y_grid, B_z(x_grid, y_grid, z_min), levels = 10)
# cl = ax.clabel(ct, inline=True, fontsize=8)
# ct = ax.contourf(x_grid, y_grid, B_z(x_grid, y_grid, z_min), levels = 10)

y_grid, z_grid = np.meshgrid(y, z)
# ct = ax.contourf(y_grid, z_grid, B_x(x_min, y_grid, z_grid), levels = 50)
im = ax.imshow(B_x(x_min, y_grid, z_grid), \
               extent = [y_min, y_max, z_min, z_max], \
               origin = 'lower')
cb = fig.colorbar(im)

# n_seeds = 21
# start_points =  np.zeros((n_seeds, 2))
# start_points[:, 0] = np.linspace(y_min, y_max, n_seeds)
# y_grid, z_grid = np.meshgrid(y, z)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ct = ax.streamplot(y_grid, z_grid, B_y(0, y_grid, z_grid), B_z(0, y_grid, z_grid), \
#                    start_points = start_points)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(mu, sol.sol(mu)[0])

plt.show()
