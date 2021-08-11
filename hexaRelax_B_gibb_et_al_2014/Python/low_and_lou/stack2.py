import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def velocity(t, coord):
    # velocity gives the velocity of the particle
    # coord gives the coordinate of the particle where
    # coord[0] = x, coord[1] = y and coord[2] = z
    return [vx_interp([coord[0], coord[1], coord[2]]), \
            vy_interp([coord[0], coord[1], coord[2]]), \
            vz_interp([coord[0], coord[1], coord[2]])]

n_points = 100
lim = 5

# Create an example of a velocity array
x = np.linspace(0, lim, n_points)
y = np.linspace(0, lim, n_points)
z = np.linspace(0, lim, n_points)
X, Y, Z = np.meshgrid(x, y, z)
vx = X
vy = Y
vz = -2 * Z

vx_interp = RegularGridInterpolator((x, y, z), vx)
vy_interp = RegularGridInterpolator((x, y, z), vy)
vz_interp = RegularGridInterpolator((x, y, z), vz)

coord0 = [2, 2, 5]
t_range = np.linspace(0, 0.2 * lim, n_points)
# sol = solve_ivp(velocity, [0, 2 * lim], t_eval = t_range, y0 = coord0)
