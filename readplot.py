# -*- coding: utf-8 -*-
#! /bin/env python3
# last updated : 2023/12/29 18:58:12
#
# Example of reading and ploting data
#
# 


import numpy as np

def read( fn ):
    """
    Read uniform grid data
    """
    from scipy.io import FortranFile
    NGH=2                       # size of ghost cells

    f = FortranFile(fn, 'r')

    # read header data
    dt = np.dtype([ ("time", "<d"), ("step", "<i")])
    chunk = f.read_record(dtype=dt)
    time = chunk[0]["time"]
    step = chunk[0]["step"]

    dt = np.dtype([("imax", "<i"), ("jmax", "<i"), ("kmax", "<i"), ("mmax", "<i")])
    chunk = f.read_record(dtype=dt)
    imax = chunk[0]["imax"]
    jmax = chunk[0]["jmax"]
    kmax = chunk[0]["kmax"]
    mmax = chunk[0]["mmax"]

    # size of arrays
    ni = imax+1+2*NGH
    nj = jmax+1+2*NGH
    nk = kmax+1+2*NGH
    nm = mmax + 1

    # example for reading multiple arrays in one record, then remove ghost cells.
    dt = np.dtype([("x", "f8", ni), ("y", "f8", nj), ("z", "f8", nk)])
    chunk = f.read_reals(dtype=dt)
    x = chunk[0]["x"][NGH:-NGH]
    y = chunk[0]["y"][NGH:-NGH]
    z = chunk[0]["z"][NGH:-NGH]

    # example for reading one array in one record, then reforming and removing ghost cells
    v = f.read_reals(dtype=np.float64).reshape((ni, nj, nk, nm),order="F")[NGH:-NGH,NGH:-NGH,NGH:-NGH,:]

    f.close()
    return x, y, z, v, time, step

def plot_2d_cont(x, y, v2d):
    """
    Plot 2D contour, v2d, in the x-y plane.
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    cs = ax.contour(x, y, v2d.T, levels=30) # note: transpose v2d
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show(block=False)
    
def plot_2d_pix(x, y, v2d):
    """
    Plot 2D pixel map, v2d, in the x-y plane.
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    extent = [x.min()-dx/2, x.max()+dx/2, y.min()-dy/2, y.max()+dy/2]
    im = ax.imshow(v2d.T, origin='lower', extent=extent) # note: transpose v2d
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show(block=False)
    
def plot_1d(x, v1d):
    """
    Plot 1D data, v1d as a function of x
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    line, = ax.plot(x, v1d)
    ax.set_xlabel('x')
    ax.set_ylabel('v')
    plt.show(block=False)
    

if __name__ == "__main__":
    """
    An example:
    If you use ipython, just type this command.
    %run readplot
    """
    import os.path
    import glob
    dir='DATA'
    prefix = 'st'
    suffix = '.d'
    fn = glob.glob(os.path.join(dir, prefix+'*'+suffix))[-1]
    print(fn)
    x, y, z, v, time, step = read(fn)
    
    if x.size != 1 and y.size != 1 and z.size == 1:
        v2d = v[:,:, 0, 0]           # 2D slice
        plot_2d_cont(x, y, v2d)
        plot_2d_pix(x, y, v2d)

    elif x.size != 1 and y.size == 1 and z.size == 1:
        v1d = v[:,0, 0, 0]           # 1D slice
        plot_1d(x, v1d)
