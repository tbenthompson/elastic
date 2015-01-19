import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt

def easy_plot_x(vertices, data):
    plt.plot(vertices[:, 0], data, 'b.-')

def quiver_plot(vertices, datax, datay):
    skip = 32
    x = vertices[::skip, 0]
    y = vertices[::skip, 1]
    opts = dict(
        scale = 20.0,
        minshaft = 1,
        headwidth = 2,
        headlength = 3,
        headaxislength = 4,
        width = 0.003
    )
    plt.quiver(x, y, datax[::skip], datay[::skip], **opts)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
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

def main(filex, filey):
    fx = h5py.File(filex)
    fy = h5py.File(filey)
    facets = fx['locations']
    datax = fx['values'][:, 0]
    datay = fy['values'][:, 0]
    plot2d(facets, datax, datay)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("FAILED")
        sys.exit()
    filenamex = sys.argv[1]
    filenamey = sys.argv[2]
    main(filenamex, filenamey)
