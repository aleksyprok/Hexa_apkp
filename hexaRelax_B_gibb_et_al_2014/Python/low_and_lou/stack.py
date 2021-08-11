import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def velocity(t, coord):
    # velocity gives the velocity of the particle
    # coord gives the coordinate of the particle where
    # coord[0] = x, coord[1] = y and coord[2] = z
    return [coord[0], coord[1], -2 * coord[2]]

n_points = 100
lim = 5

coord0 = [2, 2, 5]
t = np.linspace(0, lim, n_points)
sol = solve_ivp(velocity, [0, lim], t_eval = t, y0 = coord0)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ln = ax.plot3D(sol.y[0], sol.y[1], sol.y[2])
ax.set_xlim((0, lim))
ax.set_ylim((0, lim))
ax.set_zlim((0, lim))
plt.savefig('particle_trajectory.png')
