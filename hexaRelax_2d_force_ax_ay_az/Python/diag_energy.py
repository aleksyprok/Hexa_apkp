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
mean_divb = np.array([])
max_divb_norm = np.array([])
mean_divb_norm = np.array([])
divb_flux = np.array([])
divb_diss = np.array([])
iter_no = np.array([])
min_sin_theta = np.array([])
max_sin_theta = np.array([])
mean_sin_theta = np.array([])
sigma_j = np.array([])
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
    max_divb_norm  = np.append(max_divb_norm,  line[9])
    mean_divb_norm  = np.append(mean_divb_norm,  line[10])
    divb_flux  = np.append(divb_flux,  line[11])
    divb_diss  = np.append(divb_diss,  line[12])
    iter_no  = np.append(iter_no,  line[13])
    min_sin_theta  = np.append(min_sin_theta,  line[14])
    max_sin_theta  = np.append(max_sin_theta,  line[15])
    mean_sin_theta  = np.append(mean_sin_theta,  line[16])
    sigma_j  = np.append(sigma_j,  line[17])


print("len(t) =", len(t))

poy_int = integrate.cumtrapz(sz, t, initial = 0)
q_int = integrate.cumtrapz(q, t, initial = 0)
divb_flux_int = integrate.cumtrapz(divb_flux, t, initial = 0)
divb_diss_int = integrate.cumtrapz(divb_diss, t, initial = 0)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)

ax = fig.add_subplot(231)
ax.plot(iter_no, e - e[0])
ax.set_title(r'e - $e_0$, $e_0$ = ' + '{:1.2e}'.format(e[0]))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(232)
ax.plot(iter_no, poy_int)
ax.set_title(r'$\int sz\ dt$')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(233)
ax.plot(iter_no, q_int)
ax.set_title(r'$\int q\ dt$')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(234)
ax.plot(iter_no, divb_flux_int)
ax.set_title(r"$\int$ divB flux $dt$")
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(235)
ax.plot(iter_no, divb_diss_int)
ax.set_title(r'$\int$ divB dissipation $dt$')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(336)
ax.plot(iter_no, e - e[0] + q_int + divb_diss_int)
ax.set_title(r'$e - e_0$ + $\int$ $q$ + divB dissipation $dt$')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(337)
ax.plot(iter_no, poy_int - divb_flux_int)
ax.set_title(r'$\int$ $S_z$ - divB flux $dt$')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(338)
ax.plot(iter_no, e - e[0] + q_int + divb_diss_int - (poy_int - divb_flux_int))
ax.set_title(r'Energy error')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax = fig.add_subplot(339)
ax.plot(iter_no, t)
ax.set_title(r't')
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.show()
