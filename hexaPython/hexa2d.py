import os
import numpy as np
import matplotlib.pyplot as plt

os.makedirs('data', exist_ok = True)

def ramp_up(t):
    return np.sin(np.pi * t / 2) ** 2 * (t <= 1) + 1 * (t >  1)

def bbx_driv(x, t):
    return np.cos(np.pi * x) * ramp_up(t) * np.sin(np.pi * t)

def step():
    # This subrotuine performs 1 timestep of the algorithm
    # and updates bbx, bby, bbz.
    global t, dt

    # Calculate current at cell edges
    ccx = -(bby[1:, :] - bby[:-1, :]) / dz
    ccy =  (bbx[1:, :] - bbx[:-1, :]) / dz \
        -  (bbz[:, 1:] - bbz[:, :-1]) / dx
    ccz =  (bby[:, 1:] - bby[:, :-1]) / dx

    # Calculate current at cell corners
    cx = 0.5 * (ccx[:, :-1] + ccx[:, 1:])
    cy = ccy
    cz = 0.5 * (ccz[:-1, :] + ccz[1:, :])

    # Calculate magnetic field at cell corners
    bx = 0.5  * (bbx[:-1, :  ] + bbx[1:, :  ])
    by = 0.25 * (bby[:-1, :-1] + bby[:-1, 1:] + \
                 bby[1: , :-1] + bby[1: , 1:])
    bz = 0.5  * (bbz[:  , :-1] + bbz[:  , 1:])
    B2 = bx * bx + by * by + bz * bz

    # Calculate velocity at cell corners
    vx = frc / B2 * (cy * bz - cz * by)
    vy = frc / B2 * (cz * bx - cx * bz)
    vz = frc / B2 * (cx * by - cy * bx)

    # Calculate dt
    modb = np.sqrt(B2)
    dt = np.amin([dx / np.amax([np.amax(vx), 1e-6]), \
                  dz / np.amax([np.amax(vz), 1e-6]), \
                  dx / np.amax(modb), \
                  dz / np.amax(modb)])
    t += dt

    # Calculate E = -v x B at cell corners
    ex = vz * by - vy * bz
    ey = vx * bz - vz * bx
    ez = vy * bx - vx * by

    # Calculate E at cell edges
    eex = 0.5 * (ex[:, :-1] + ex[:, 1:])
    eey = ey
    eez = 0.5 * (ez[:-1, :] + ez[1:, :])

    # Update bbx, bby, bbz using dB/dt = -curl(E)
    bbx[1:-1, :   ] += dt *  (eey[1:, :] - eey[:-1, :]) / dz
    bby[1:-1, 1:-1] += dt * ((eez[:, 1:] - eez[:, :-1]) / dx - \
                             (eex[1:, :] - eex[:-1, :]) / dz)
    bbz[:   , 1:-1] -= dt *  (eey[:, 1:] - eey[:, :-1]) / dx

def boundary_conditions():
    # This subrotuine imposes the boundary conditions

    # x = x_min
    bby[:, 0] = bby[:, 1]
    bbz[:, 0] = bbz[:, 1]

    # x = x_max
    bby[:, -1] = bby[:, -2]
    bbz[:, -1] = bbz[:, -2]

    # z = z_min
    bbx[0, :] = 2 * bbx_driv(xb, t) - bbx[1, :]
    bby[0, :] = 2 * bby_driv(xc, t) - bby[1, :]

    # z = z_max
    bbx[-1, :] = bbx[-2, :]
    bby[-1, :] = bby[-2, :]

def overwrite():
    # This subroutine overwrites bz at z = z_min
    bbz[0, :] = np.ones(nx + 2)

frc = 1e-3 # Frictional coefficent

nx = 128
nz = 128
nt = 1024
t = 0
t_array = np.zeros(1)
dt_snapshots = 0.1
snapshot_no = 0

# Create grid
x_min = -0.5
x_max =  0.5
dx = (x_max - x_min) / nx
xb = np.linspace(x_min, x_max, nx + 1)
xc = np.linspace(x_min - 0.5 * dx, x_max + 0.5 * dx, nx + 2)

z_min = 0
z_max = 1
dz = (z_max - z_min) / nz
zb = np.linspace(z_min, z_max, nz + 1)
zc = np.linspace(z_min - 0.5 * dz, z_max + 0.5 * dz, nz + 2)

# Initial condition
bbx = np.zeros((nz + 2, nx + 1))
bby = np.zeros((nz + 2, nx + 2))
bbz = np.ones( (nz + 1, nx + 2))

# Output inital arrays
np.save('data/bbx' + "{:04d}".format(snapshot_no),  bbx)
np.save('data/bby' + "{:04d}".format(snapshot_no),  bby)
np.save('data/bbz' + "{:04d}".format(snapshot_no),  bbz)
np.save('data/t' + "{:04d}".format(snapshot_no),  t)

for n in range(nt):
    step()
    boundary_conditions()
    overwrite()
    if t > t_array[snapshot_no] + dt_snapshots:
        t_array = np.append(t_array, t)
        snapshot_no += 1
        np.save('data/bbx' + "{:04d}".format(snapshot_no),  bbx)
        np.save('data/bby' + "{:04d}".format(snapshot_no),  bby)
        np.save('data/bbz' + "{:04d}".format(snapshot_no),  bbz)
        np.save('data/t'   + "{:04d}".format(snapshot_no),  t)
    # np.save('data/bbx' + "{:04d}".format(n + 1),  bbx)
    # np.save('data/bby' + "{:04d}".format(n + 1),  bby)
    # np.save('data/bbz' + "{:04d}".format(n + 1),  bbz)
    # np.save('data/t' + "{:04d}".format(n + 1),  t)
