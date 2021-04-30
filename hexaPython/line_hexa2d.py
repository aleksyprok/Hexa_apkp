import os
import numpy as np
import matplotlib.pyplot as plt

# Make directories
os.makedirs('figures/line/bbx', exist_ok=True)
os.makedirs('figures/line/bby', exist_ok=True)
os.makedirs('figures/line/bbz', exist_ok=True)

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

for n in range(num_files):

    bbx = np.load('data/bbx' + "{:04d}".format(n) + '.npy')
    bby = np.load('data/bby' + "{:04d}".format(n) + '.npy')
    bbz = np.load('data/bbz' + "{:04d}".format(n) + '.npy')
    t  = np.load('data/t'  + "{:04d}".format(n) + '.npy')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(zc, bbx[:, nx // 2])
    ax.set_xlabel('z')
    ax.set_title(r'$B_x$ at $t =$' + '{:10e}'.format(t))
    fig.savefig('figures/line/bbx/' + "{:04d}".format(n) + '.png', bbox_inches = 'tight')
    plt.close(fig)
    #
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(zc, bby[:, nx // 2])
    ax.set_xlabel('z')
    ax.set_title(r'$B_y$ at $t =$' + '{:10e}'.format(t))
    fig.savefig('figures/line/bby/' + "{:04d}".format(n) + '.png', bbox_inches = 'tight')
    plt.close(fig)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(zb, bbz[:, nx // 2])
    ax.set_xlabel('z')
    ax.set_title(r'$B_z$ at $t =$' + '{:10e}'.format(t))
    fig.savefig('figures/line/bbz/' + "{:04d}".format(n) + '.png', bbox_inches = 'tight')
    plt.close(fig)
