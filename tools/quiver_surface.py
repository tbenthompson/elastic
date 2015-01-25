import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt

def easy_plot_x(vertices, data):
    plt.plot(vertices[:, 0], data, 'b.-')

def quiver_plot(vertices, datax, datay):
    skip = vertices.shape[0] / 64
    x = vertices[::skip, 0]
    y = vertices[::skip, 1]
    scale = 12.0 * (np.mean(np.abs(datax)) + np.mean(np.abs(datay)))
    opts = dict(
        scale = scale,
        minshaft = 1,
        headwidth = 2,
        headlength = 3,
        headaxislength = 4,
        width = 0.003
    )
    plt.quiver(x, y, datax[::skip], datay[::skip], **opts)

    view_factor = 0.2
    plt.xlim(np.min(x) + (np.min(x) - np.mean(x)) * view_factor,
              np.max(x) + (np.max(x) - np.mean(x)) * view_factor)
    plt.ylim(np.min(y) + (np.min(y) - np.mean(y)) * view_factor,
              np.max(y) + (np.max(y) - np.mean(y)) * view_factor)
    plt.show()

def plot2d(facets, datax, datay):
    x_index = [0, 2]
    y_index = [1, 3]
    vertices = np.array([
        facets[:, x_index].flatten(),
        facets[:, y_index].flatten()
    ]).T

    x = vertices[:, 1]
    # easy_plot_x(vertices, data[:, values_dim])
    quiver_plot(vertices, datax, datay)

    plt.show()

def main(filename):
    f = h5py.File(filename)
    facets = f['locations']
    datax = f['values0'][:,0]
    datay = f['values1'][:,0]
    plot2d(facets, datax, datay)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("FAILED")
        sys.exit()
    main(sys.argv[1])
