import numpy as np
from scipy.linalg import solve_banded

# !============================================================
# ! Solutions to a system of tridiagonal linear equations C*x=b
# ! Method: the Thomas method
# ! Alex G. November 2009
# !-----------------------------------------------------------
# ! input ...
# ! c(n,3) - array of coefficients for matrix C
# ! b(n)   - vector of the right hand coefficients b
# ! n      - number of equations
# ! output ...
# ! x(n)   - solutions
# ! comments ...
# ! the original arrays c(n,3) and b(n) will be destroyed
# ! during the calculation
# !===========================================================

nx = 16

c = np.zeros((3, nx))
b = np.zeros(nx)
c[0, :] = 0.5
c[1, :] = 1.0
c[2, :] = -0.5
b[-1] = 10

x = np.zeros(nx)

# step 1: forward elimination
for ix in range(1, nx):
    coeff = c[0, ix] / c[1, ix-1]
    c[1, ix] = c[1, ix] - coeff * c[2, ix-1]
    b[ix] = b[ix] - coeff * b[ix-1]

# step 2: back substitution
x[nx-1] = b[nx-1] / c[1, nx-1]
for ix in range(nx-2, -1, -1):
    x[ix] = (b[ix]- c[2, ix] * x[ix+1]) / c[1, ix]

print(x)

c = np.zeros((3, nx))
b = np.zeros(nx)
c[0, :] = -0.5
c[1, :] = 1.0
c[2, :] = 0.5
b[-1] = 10

print(solve_banded((1,1), c, b))
