import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

t = np.array([])
dt = np.array([])
e  = np.array([])
sz = np.array([])
q = np.array([])
err = np.array([])
min_divb = np.array([])
max_divb = np.array([])
mean_divb = np.array([])
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
    err  = np.append(err,  line[5])
    min_divb  = np.append(min_divb,  line[6])
    max_divb  = np.append(max_divb,  line[7])
    mean_divb  = np.append(mean_divb,  line[8])

print("len(t) =", len(t))

base = integrate.cumtrapz(sz, t, initial = 0)
diss = integrate.cumtrapz(q, t, initial = 0)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(331)
ax.plot(t, e - e[0])
ax.set_title(r'e - $e_0$, $e_0$ = ' + '{:1.2e}'.format(e[0]))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(332)
ax.plot(t, base)
ax.set_title(r'$\int sz\ dt$')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(333)
ax.plot(t, diss)
ax.set_title(r'$\int q\ dt$')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(334)
ax.plot(t, e - e[0] + diss)
ax.set_title(r'$e - e_0 + \int q\ dt$')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(335)
ax.plot(t, e - e[0] + diss - base)
ax.set_title(r'$e - e_0 + \int q\ dt - \int sz\ dt$')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(336)
ax.plot(t)
ax.set_title('t')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(337)
ax.plot(t, err)
ax.set_title(r'$\int |\mathbf{B}_{num} - \mathbf{B}_{ana}|^2 dV$')
lowerlim = 0
upperlim = ax.get_ylim()[1]
ax.set_ylim(lowerlim, upperlim)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(338)
ax.plot(t, max_divb)
ax.set_title('Max divB')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(339)
ax.plot(t, mean_divb)
ax.set_title('Mean divB')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.show()
