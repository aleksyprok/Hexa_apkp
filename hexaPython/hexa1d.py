import numpy as np
import matplotlib.pyplot as plt

def Bx_driv(t):
    return np.sin(t)

def step():
    # This subrotuine performs 1 timestep of the algorithm
    # and upates Bx and By.
    global t

    # Calculate current
    jx = -(By[1:] - By[:-1]) / dz
    jy =  (Bx[1:] - Bx[:-1]) / dz

    # Calculate velocity
    Bx_avg = 0.5 * (Bx[:-1] + Bx[1:])
    By_avg = 0.5 * (By[:-1] + By[1:])
    B2   = Bx_avg * Bx_avg + By_avg * By_avg + Bz * Bz
    vx =  frc * jy * Bz / B2
    vy = -frc * jx * Bz / B2
    vz =  frc * (jx * By_avg - jy * Bx_avg) / B2

    # Calculate dt
    modb = np.sqrt(B2 + Bz ** 2)
    # dt = np.amin([dz / np.amax([np.amax(vz), 1e-6]), dz / np.amax(modb)])
    dt = 1e-2
    t += dt
    t_array[n+1] = t

    # Calculate v x B
    cx = vy * Bz     - vz * By_avg
    cy = vz * Bx_avg - vx * Bz

    # Calculate curl of v x B
    ccx = -(cy[1:] - cy[:-1]) / dz
    ccy =  (cx[1:] - cx[:-1]) / dz

    # Update Bx, By
    Bx[1:-1] += dt * ccx
    By[1:-1] += dt * ccy

def boundary_conditions1():
    Bx[0] = -Bx[1] + 2 * Bx_driv(t)
    By[0] = By[1]
    Bx[-1] = Bx[-2]
    By[-1] = By[-2]

def boundary_conditions2():
    Bx[0] = 1e6
    Bx[1] = Bx_array1[n + 1, 1]
    By[0] = By[1]
    Bx[-1] = Bx[-2]
    By[-1] = By[-2]

Bz = 1
frc = 1e-2 # Frictional coefficent

nz = 64
nt = 1024
t_array = np.zeros(nt + 1)

# Create grid
z_min = 0
z_max = 1
dz = (z_max - z_min) / nz
# zb gives z at the cell boundaries
zb = np.linspace(z_min, z_max, nz + 1)
# zc gives z at the cell centres,
# where z[0] and z[-1] gives the cell centres of the ghost cells
zc = np.linspace(z_min - 0.5 * dz, z_max + 0.5 * dz, nz + 2)

# Initialise data output
Bx_array = np.zeros((nt + 1, nz + 2))
By_array = np.zeros((nt + 1, nz + 2))
Bx_array1 = np.zeros((nt + 1, nz + 2))
Bx_array2 = np.zeros((nt + 1, nz + 2))

# Run simulation with boundary condition 1
# Initial condition
Bx = np.zeros(nz + 2)
By = np.zeros(nz + 2)
t = 0
for n in range(nt):
    step()
    boundary_conditions1()
    # Output data
    Bx_array[n + 1, :] = Bx
    By_array[n + 1, :] = By
Bx_array1[:, :] = Bx_array[:, :]

# Run simulation with boundary condition 2
# Initial condition
Bx = np.zeros(nz + 2)
By = np.zeros(nz + 2)
t = 0
for n in range(nt):
    step()
    boundary_conditions2()
    # Output data
    Bx_array[n + 1, :] = Bx
    By_array[n + 1, :] = By
Bx_array2[:, :] = Bx_array[:,:]

# Plot Bx along z movie
for n in range(nt + 1):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(zc[0:10], Bx_array1[n, 0:10], '-+', label = 'Boundary condition 1')
    ax.plot(zc[0:10], Bx_array2[n, 0:10], '-+', label = 'Boundary condition 2')
    ax.set_xlabel('z')
    ax.set_title(r'$B_x$ at $t =$' + '{:06.2f}'.format(t_array[n]))
    ax.set_ylim(-1.5, 1.5)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, bbox_to_anchor = (0.7, 1.25))
    fig.savefig('figures/Bx' + "{:04d}".format(n) + '.png', bbox_inches = 'tight')
    plt.close(fig)
