import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import glob

nx = 128
ny = 128
nz = 128
dx = np.float64(6 / nx)
dy = np.float64(6 / ny)
dz = np.float64(6 / nz)

xc = np.arange(dx / 2, 6, dx)
yc = np.arange(dy / 2, 6, dy)
zc = np.arange(dz / 2, 6, dz)
xb = np.arange(0, 6 + dx /2, dx)
yb = np.arange(0, 6 + dy /2, dy)
zb = np.arange(0, 6 + dz /2, dz)

aax = np.zeros((nz+1, ny+1, nx), dtype = 'float32')
aay = np.zeros((nz+1, ny, nx+1), dtype = 'float32')
aaz = np.zeros((nz, ny+1, nx+1), dtype = 'float32')

for ix in range(nx):
    for iy in range(ny+1):
        for iz in range(nz+1):
            aax[iz, iy, ix] = 0

for ix in range(nx+1):
    for iy in range(ny):
        for iz in range(nz+1):
            aay[iz, iy, ix] = xb[ix]

for ix in range(nx+1):
    for iy in range(ny+1):
        for iz in range(nz):
            aaz[iz, iy, ix] = 0

f = FortranFile('setup_files/simple_ic', 'w')
f.write_record(np.int32(1))
f.write_record(aax)
f.write_record(aay)
f.write_record(aaz)
