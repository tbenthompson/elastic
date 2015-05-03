import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import griddata

def main(filename, which_field, scale):
    f = h5py.File(filename)
    locs = f['locations']
    x = locs[:, 0]
    y = locs[:, 1]
    data = f['values' + str(which_field)][:, 0]
    if scale == 1:
        data = np.log10(np.abs(data))
    elif scale == 2:
        data = np.sign(data) * (np.abs(data) ** (1.0 / 2.0))

    xrange = [np.min(x), np.max(x)]
    yrange = [np.min(y), np.max(y)]

    view_factor = 0.2
    xlims = [xrange[0] + (np.min(x) - np.mean(x)) * view_factor,
             xrange[1] + (np.max(x) - np.mean(x)) * view_factor]
    ylims = [yrange[0] + (np.min(y) - np.mean(y)) * view_factor,
             yrange[1] + (np.max(y) - np.mean(y)) * view_factor]

    nx = 500
    ny = 500

    xi = np.linspace(xrange[0] + 0.0001, xrange[1] - 0.0001, nx)
    yi = np.linspace(yrange[0] + 0.0001, yrange[1] - 0.0001, ny)
    # grid the data.
    zi = griddata((x, y), data, (xi[None,:], yi[:,None]))

    data_min = np.min(data)
    data_max = np.max(data)
    # levels = np.linspace(-1.0, 0.5, 20)
    levels = np.linspace(data_min, data_max, 20)
    plt.figure()
    plt.contourf(xi, yi, zi, levels = levels, extend = 'both')
    plt.colorbar()
    # plt.contour(xi, yi, zi, linestyles = 'solid', colors = '#000000', levels = levels, extend = 'both')
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    title_str = filename + ' component: ' + str(which_field)
    if scale == 1:
        title_str += ' LOG SCALE'
    elif scale == 2:
        title_str += ' CUBE ROOT SCALE'
    plt.title(title_str)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
    plt.show()
