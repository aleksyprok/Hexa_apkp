import numpy as np
from scipy.integrate import solve_bvp
from scipy.misc import derivative
import matplotlib.pyplot as plt
from scipy.io import FortranFile

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
    return -sol.sol(mu)[1] / r ** 3

def B_theta(r, theta):
    mu = np.cos(theta)
    return sol.sol(mu)[0] / (r ** 3 * np.sin(theta))

def B_phi(r, theta):
    mu = np.cos(theta)
    return sol.p[0] * sol.sol(mu)[0] ** 2 / (r ** 3 * np.sin(theta))

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

lx = 0.0
lz = -0.25
Phi = 4 * np.pi / 10

nx = 128
x_min = -0.8
x_max = 1.2
dx = (x_max - x_min) / nx
x = np.linspace(x_min, x_max, nx)
xb = np.linspace(x_min, x_max, nx + 1)
xc = np.linspace(x_min - dx / 2, x_max + dx / 2, nx + 2)

ny = 128
y_min = -1
y_max = 1
dy = (y_max - y_min) / ny
y = np.linspace(y_min, y_max, ny)
yb = np.linspace(y_min, y_max, ny + 1)
yc = np.linspace(y_min - dy / 2, y_max + dy / 2, ny + 2)

nz = 128
z_min = 1
z_max = 3
dz = (z_max - z_min) / nz
z = np.linspace(z_min, z_max, nz)
zb = np.linspace(z_min, z_max, nz + 1)
zc = np.linspace(z_min - dz / 2, z_max + dz / 2, nz + 2)

X, Y, Z = np.meshgrid(xc, yc, zb, indexing = 'ij')
bbz = B_z(X, Y, Z)
b_norm = 500 / np.max(np.abs(bbz))
bbz = B_z(X, Y, Z) * b_norm

X, Y, Z = np.meshgrid(xb, yc, zc, indexing = 'ij')
bbx = B_x(X, Y, Z) * b_norm

X, Y, Z = np.meshgrid(xc, yb, zc, indexing = 'ij')
bby = B_y(X, Y, Z) * b_norm

f = FortranFile('setup_files/low_and_low_nlff.dat', 'w')
f.write_record(bbx.T)
f.write_record(bby.T)
f.write_record(bbz.T)
f.close()
print(bbx[:, ny//2, nz//2])
