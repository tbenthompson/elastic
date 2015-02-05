import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import sys

def main(filename):
    f = h5py.File(filename)
    locs = f['locations']
    x = locs[:, 0]
    y = locs[:, 1]
    datax = f['values0'][:,0]
    datay = f['values1'][:,0]

    scale = 50.0 * (np.mean(np.abs(datax)) + np.mean(np.abs(datay)))
    opts = dict(
        scale = scale,
        minshaft = 1,
        headwidth = 2,
        headlength = 3,
        headaxislength = 4,
        width = 0.003
    )

    print datax[x < -0.99]
    print datay[x < -0.99]


    plt.quiver(x, y, datax, datay, **opts)
    view_factor = 0.2
    plt.xlim(np.min(x) + (np.min(x) - np.mean(x)) * view_factor,
              np.max(x) + (np.max(x) - np.mean(x)) * view_factor)
    plt.ylim(np.min(y) + (np.min(y) - np.mean(y)) * view_factor,
              np.max(y) + (np.max(y) - np.mean(y)) * view_factor)
    plt.title(filename)
    plt.show()

    # define grid.
    # xi = np.linspace(-1.0,1.0,100)
    # yi = np.linspace(-1.0,1.0,100)
    # # grid the data.
    # from scipy.interpolate import griddata
    # zi = griddata((x, y), datay, (xi[None,:], yi[:,None]), method='cubic')
    # plt.contourf(xi, yi, zi)
    # plt.contour(xi, yi, zi)

if __name__ == "__main__":
    main(sys.argv[1])
