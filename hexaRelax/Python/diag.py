import numpy as np
import matplotlib.pyplot as plt

t = np.array([])
e  = np.array([])
q = np.array([])
min_divB = np.array([])
max_divB = np.array([])
mean_divB = np.array([])
file = open('run1/diagnostic', 'r')
n = 0
while True:
    n += 1
    line = file.readline()
    if n % 1000 == 0: print(n)
    if not line:
        break
    line = np.float64(line.split())
    # line = line.split())
    t = np.append(t,  line[0])
    e = np.append(e,  line[1])
    q = np.append(q,  line[2])
    min_divB  = np.append(min_divB,  line[3])
    max_divB  = np.append(max_divB,  line[4])
    mean_divB  = np.append(mean_divB,  line[5])

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)

ax = fig.add_subplot(231)
ax.plot(e - e[0])
ax.set_title(r'Magnetic energy (ergs), $E_0$ = ' + '{:1.2e}'.format(e[0]))

ax = fig.add_subplot(232)
ax.plot(t)
ax.set_title('t')

ax = fig.add_subplot(233)
ax.plot((t[1:] - t[:-1]))
ax.set_title('dt')

ax = fig.add_subplot(234)
ax.plot(q)
ax.set_title('Frictional dissipation (ergs/s)')

ax = fig.add_subplot(235)
ax.plot(max_divB)
ax.set_title('Max divB')

ax = fig.add_subplot(236)
ax.plot(mean_divB)
ax.set_title('Mean divB')

plt.show()
