import numpy as np
from scipy.io import FortranFile

dir = '../hexaf90/run1/'

# Get nx, ny, nz
f = open(dir + 'param1')
num_hex_cells = int(f.readline())
nx_ny_nz_string = f.readline()
nx = int(nx_ny_nz_string.split()[0])
ny = int(nx_ny_nz_string.split()[1])
nz = int(nx_ny_nz_string.split()[2])

# Get nt
f = open(dir + 'run1_setup')
vsetup = f.readline()
nmajor_string = f.readline()
nmajor = int(list(filter(str.isdigit, nmajor_string))[0])
print(nmajor)

# Get aax, aay, aaz
root = 'run1_00001p'
f = FortranFile(dir + root)
opt = f.read_ints()[0]
aax = f.read_reals(dtype = 'float32').reshape((nz + 1, ny + 1, nx    ))
aay = f.read_reals(dtype = 'float32').reshape((nz + 1, ny    , nx + 1))
aaz = f.read_reals(dtype = 'float32').reshape((nz    , ny + 1, nx + 1))

# Get aax0, aay0
root = 'run1_evolve'
f = FortranFile(dir + root)
opt = f.read_ints()[0]
aax0 = np.zeros((nmajor, ny + 1, nx    ), dtype = 'float32')
aay0 = np.zeros((nmajor, ny    , nx + 1), dtype = 'float32')
for n in range(nmajor):
    aax0[n,:,:] = f.read_reals(dtype = 'float32').reshape((ny + 1, nx    ))
    aay0[n,:,:] = f.read_reals(dtype = 'float32').reshape((ny    , nx + 1))

print(aax0[1,0:2,0:2])
print( aax[0,0:2,0:2])
print(' ')
print(aay0[1,0:2,0:2])
print( aay[0,0:2,0:2])
