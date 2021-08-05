import numpy as np
from streamtracer import StreamTracer, VectorGrid

nsteps = 10000
step_size = 0.1
tracer = StreamTracer(nsteps, step_size)

field = np.ones((10, 10, 10, 3))
grid_spacing = [1, 2, 1]
grid = VectorGrid(field, grid_spacing)

seeds = np.array([[0, 0, 0], [0, 0, 1]])
tracer.trace(seeds, grid)
