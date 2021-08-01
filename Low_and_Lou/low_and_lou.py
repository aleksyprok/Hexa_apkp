import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

def dY_dmu(mu, Y, p):
    # Let Q = (d P) / (d mu)
    # Y[0] = P
    # Y[1] = Q
    n = 1
    a = p[0]
    dP_dmu = Y[1]
    dQ_dmu = n * (n + 1) * Y[0] + a * a * (1 + n) / n * Y[0] ** (1 + 2 / n)
    dQ_dmu = -dQ_dmu / (1 - mu * mu)
    return np.array([dP_dmu, dQ_dmu])

def bcs(Y_a, Y_b, p):
	return np.array([Y_a[0], Y_b[0], Y_a[1] - 10])

mu_min = -1 + 1e-12
mu_max = 1 - 1e-12
mu_temp = np.linspace(mu_min, mu_max, 5)
Y = np.zeros((2, mu_temp.size))
Y[0, 1] = 1
Y[0, 3] = -1
# mu_temp = np.linspace(mu_min, mu_max, 7)
# Y = np.zeros((2, mu_temp.size))
# Y[0, 1] = 1
# Y[0, 3] = -1
# Y[0, 5] = 1

n_mu = 1024
mu = np.linspace(mu_min, mu_max, n_mu)

sol = solve_bvp(dY_dmu, bcs, mu_temp, Y, p = [np.sqrt(0.425)], tol = 1e-10)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(mu, sol.sol(mu)[0])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(np.arccos(mu), sol.sol(mu)[0])

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(mu, sol.sol(mu)[1])

print(sol.p[0] ** 2)