import os
import numpy as np
import matplotlib.pyplot as plt

# Make directories
os.makedirs('figures/surface/bbx', exist_ok=True)
os.makedirs('figures/surface/bby', exist_ok=True)
os.makedirs('figures/surface/bbz', exist_ok=True)
os.makedirs('figures/surface/divB', exist_ok=True)

# Count number of files
num_files = 0
filenames = os.listdir('data')
for filename in filenames:
    if filename.startswith('bbx'):
        num_files += 1

# Get nx and nz
bbx = np.load('data/bbx' + "{:04d}".format(0) + '.npy')
nx = len(bbx[0, :]) - 1
nz = len(bbx[:, 0]) - 2

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

snapshot = num_files // 2

bbx = np.load('data/bbx' + "{:04d}".format(snapshot) + '.npy')
bby = np.load('data/bby' + "{:04d}".format(snapshot) + '.npy')
bbz = np.load('data/bbz' + "{:04d}".format(snapshot) + '.npy')

# X, Z = np.meshgrid(xb, zc)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Z, bbx)
#
# X, Z = np.meshgrid(xc, zc)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Z, bby)

# X, Z = np.meshgrid(xc, zb)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Z, bbz)

for n in range(num_files):

    bbx = np.load('data/bbx' + "{:04d}".format(n) + '.npy')
    bby = np.load('data/bby' + "{:04d}".format(n) + '.npy')
    bbz = np.load('data/bbz' + "{:04d}".format(n) + '.npy')
    t  = np.load('data/t'  + "{:04d}".format(n) + '.npy')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Z = np.meshgrid(xb, zc)
    ax.plot_surface(X, Z, bbx)
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title(r'$B_x$ at $t =$' + '{:10e}'.format(t))
    fig.savefig('figures/surface/bbx/' + "{:04d}".format(n) + '.png', bbox_inches = 'tight')
    plt.close(fig)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Z = np.meshgrid(xc, zc)
    ax.plot_surface(X, Z, bby)
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title(r'$B_y$ at $t =$' + '{:10e}'.format(t))
    fig.savefig('figures/surface/bby/' + "{:04d}".format(n) + '.png', bbox_inches = 'tight')
    plt.close(fig)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Z = np.meshgrid(xc, zb)
    ax.plot_surface(X, Z, bbz)
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title(r'$B_z$ at $t =$' + '{:10e}'.format(t))
    fig.savefig('figures/surface/bbz/' + "{:04d}".format(n) + '.png', bbox_inches = 'tight')
    plt.close(fig)

    divB = (bbx[1:-1, 1:] - bbx[1:-1, :-1]) / dx \
         + (bbz[1:, 1:-1] - bbz[:-1, 1:-1]) / dz
    print('t = ' + '{:5e}'.format(t) + ', max divB = ' + '{:5e}'.format(np.amax(divB[1:,:])))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Z = np.meshgrid(xc[1:-1], zc[1:-1])
    ax.plot_surface(X, Z, divB)
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title(r'$\mathbf{\nabla \cdot B}$ at $t =$' + '{:10e}'.format(t))
    fig.savefig('figures/surface/divB/' + "{:04d}".format(n) + '.png', bbox_inches = 'tight')
    plt.close(fig)
