import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

t = np.array([])
dt = np.array([])
e  = np.array([])
sz = np.array([])
q = np.array([])
min_diva = np.array([])
max_diva = np.array([])
mean_diva = np.array([])
file = open('run1/diagnostic', 'r')
n = 0
while True:
    n += 1
    line = file.readline()
    if n % 1000 == 0: print(n)
    if not line:
        break
    line = np.float64(line.split())
    t = np.append(t,  line[0])
    dt = np.append(dt,  line[1])
    e = np.append(e,  line[2])
    sz = np.append(sz,  line[3])
    q = np.append(q,  line[4])
    min_diva  = np.append(min_diva,  line[5])
    max_diva  = np.append(max_diva,  line[6])
    mean_diva  = np.append(mean_diva,  line[7])

print(t[0])

# base = np.cumsum(sz * dt)
# diss = np.cumsum(q * dt)
base = integrate.cumtrapz(sz, t, initial = 0)
diss = integrate.cumtrapz(q, t, initial = 0)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1]
fig.set_size_inches(fig_size)

ax = fig.add_subplot(131)
ax.plot(t, (e - e[0]) / e[0])
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(132)
ax.plot(t, e - e[0])
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(133)
ax.plot(t, e)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.show()
