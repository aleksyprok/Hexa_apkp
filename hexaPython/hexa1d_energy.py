import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import integrate

os.makedirs('figures', exist_ok = True)

def ramp_up(t):
    return np.sin(np.pi * t) ** 2 * (t <= 1.0) + 1.0 * (t > 1.0)

def bbx_driv(t):
    return np.sin(np.pi * t) * ramp_up(t)

def step():
    # This subrotuine performs 1 timestep of the algorithm
    # and upates bbx.
    global t, bbx, eey, bx, vx, vz

    # Calculate current
    ccy = (bbx[1:] - bbx[:-1]) / dz

    # Calculate velocity
    bx = 0.5 * (bbx[:-1] + bbx[1:])
    bb = bx * bx + bz * bz
    vx =  frc * ccy * bz / bb
    vz = -frc * ccy * bx / bb

    # Calculate t
    t += dt
    t_array[n+1] = t

    # Calculate E = - v x B
    eey = vx * bz - vz * bx

    # Overwrite E at z = 0

    # Update bbx from dBx / dt = - curl(E)
    bbx[1:-1] += dt * (eey[1:] - eey[:-1]) / dz

def boundary_conditions():
    bbx[0] = -bbx[1] + 2 * bbx_driv(t)
    bbx[-1] = -bbx[-2]

bz = 1
frc = 1e-2 # Frictional coefficent

# Create grid
nz = 256
z_min = 0
z_max = 1
dz = (z_max - z_min) / nz
zb = np.linspace(z_min, z_max, nz + 1)
zc = np.linspace(z_min - 0.5 * dz, z_max + 0.5 * dz, nz + 2)

dt = 0.05 * dz
t_max = 4
nt = int(np.ceil(t_max / dt))
t_array = np.zeros(nt + 1)

# Initialise data output
bbx_array = np.zeros((nt + 1, nz + 2))
bx_array = np.zeros((nt + 1, nz + 1))
eey_array = np.zeros((nt + 1, nz + 1))
vx_array = np.zeros((nt + 1, nz + 1))
vz_array = np.zeros((nt + 1, nz + 1))

# Initial condition
bbx = np.zeros(nz + 2)
t = 0
for n in range(nt):
    step()
    boundary_conditions()
    # Output data
    bbx_array[n + 1, :] = bbx
    bx_array[n + 1, :] = bx
    eey_array[n + 1, :] = eey
    vx_array[n + 1, :] = vx
    vz_array[n + 1, :] = vz

e_tot = np.sum(bbx_array[:, 1:-1] * bbx_array[:, 1:-1], axis = 1) / 2 * dz

poy_flux = -eey_array[:, 0] * bx_array[:, 0]
poy_tot = integrate.cumtrapz(poy_flux, t_array, initial = 0)

vx = 0.5 * (vx_array[:, :-1] + vx_array[:, 1:])
vz = 0.5 * (vz_array[:, :-1] + vz_array[:, 1:])
vv = vx * vx + vz * vz
bb = bbx_array[:, 1:-1] * bbx_array[:, 1:-1] + bz * bz
q_diss = np.sum(vv * bb, axis = 1) * dz / frc
q_tot = integrate.cumtrapz(q_diss, t_array, initial = 0)

# fig = plt.figure()
# fig_size = fig.get_size_inches()
# fig_size[0] = fig_size[0] * 3
# fig_size[1] = fig_size[1] * 2
# fig.set_size_inches(fig_size)
#
# ax = fig.add_subplot(231)
# ax.plot(t_array, e_tot)
# ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#
# ax = fig.add_subplot(232)
# ax.plot(t_array, poy_tot - q_tot)
# ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#
# ax = fig.add_subplot(233)
# ax.plot(t_array, e_tot - (poy_tot - q_tot))
# ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#
# ax = fig.add_subplot(234)
# ax.plot(t_array, poy_tot)
# ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#
# ax = fig.add_subplot(235)
# ax.plot(t_array, q_tot)
# ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#
# plt.show()

# Plot bbx along z movie
dt_snapshots = 0.1
t_snapshot = 0
for n in range(0, nt + 1):

    if (t_array[n] >= t_snapshot):
        t_snapshot += dt_snapshots

        fig = plt.figure()
        fig_size = fig.get_size_inches()
        fig_size[0] = fig_size[0] * 3
        fig_size[1] = fig_size[1] * 1
        fig.set_size_inches(fig_size)

        ax = fig.add_subplot(131)
        ax.plot(zc, bbx_array[n, :], '-+')
        ax.set_ylim(-1.1, 1.1)
        ax.set_xlabel('z')
        ax.set_title(r'$B_x$ at $t =$' + '{:07.3f}'.format(t_array[n]))

        ax = fig.add_subplot(132)
        ax.plot(zb, vx_array[n, :], '-+')
        ax.set_ylim(-0.3, 0.3)
        ax.set_xlabel('z')
        ax.set_title(r'$v_x$ at $t =$' + '{:07.3f}'.format(t_array[n]))

        ax = fig.add_subplot(133)
        ax.plot(zb, vz_array[n, :], '-+')
        ax.set_ylim(-0.3, 0.3)
        ax.set_xlabel('z')
        ax.set_title(r'$v_z$ at $t =$' + '{:07.3f}'.format(t_array[n]))

        fig.savefig('figures/' + "{:05d}".format(n) + '.png', bbox_inches = 'tight')
        plt.close(fig)
