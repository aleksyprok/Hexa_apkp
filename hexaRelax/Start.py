import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

def getdata(n):

    global bbx, bby, bbz

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    bbx = file.read_reals('float32').reshape((nz + 2, ny + 2, nx + 1), order = "C")
    bby = file.read_reals('float32').reshape((nz + 2, ny + 1, nx + 2), order = "C")
    bbz = file.read_reals('float32').reshape((nz + 1, ny + 2, nx + 2), order = "C")

def get_XY(var):

    global x, y, X, Y

    if var == 'bbx':
        x = np.arange(1, nx+2)
        y = np.arange(ny+2)
        X, Y = np.meshgrid(x, y)
    elif var == 'bby':
        x = np.arange(nx+2)
        y = np.arange(1, ny+2)
        X, Y = np.meshgrid(x, y)
    elif var == 'bbz':
        x = np.arange(nx+2)
        y = np.arange(ny+2)
        X, Y = np.meshgrid(x, y)

nx = 128
ny = 128
nz = 128
