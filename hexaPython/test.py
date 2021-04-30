import numpy as np
from scipy.io import FortranFile

nx = 2
ny = 2
nz = 3

f = FortranFile('test.dat')
aax = f.read_reals(dtype = 'float32').reshape((nz + 1, ny + 1, nx))
print(aax)
# for iz in range(nz + 1):
#     for iy in range(ny + 1):
#         print(aax[iz,iy,:])
#     print(' ')
#     print(' ')
