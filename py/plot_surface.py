import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import mayavi.mlab as mlab

def easy_plot_x(vertices, data):
    plt.plot(vertices[:, 0], data, 'b.-')

def quiver_plot(vertices, data):
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
    plt.quiver(x, y, data[::skip, 0], data[::skip, 1], **opts)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.show()

def plot2d(facets, data, values_dim):
    x_index = [0, 2]
    y_index = [1, 3]
    vertices = np.array([
        facets[:, x_index].flatten(),
        facets[:, y_index].flatten()
    ]).T

    x = vertices[:, 1]
    # easy_plot_x(vertices, data[:, values_dim])
    quiver_plot(vertices, data)

    plt.show()

def plot3d(facets, data, values_dim):
    x_index = [0, 3, 6]
    y_index = [1, 4, 7]
    z_index = [2, 5, 8]
    vertices = np.array([
        facets[:, x_index].flatten(),
        facets[:, y_index].flatten(),
        facets[:, z_index].flatten()
    ]).T

    n_v = vertices.shape[0]
    faces = np.arange(n_v).reshape((n_v / 3, 3))

    mlab.triangular_mesh(vertices[:, 0], vertices[:, 1], vertices[:, 2],
        faces[:,:], scalars = data[:, values_dim])
    mlab.show()

def main(filename, values_dim):
    f = h5py.File(filename)
    facets = f['locations']
    data = f['values']
    if (facets.shape[1] == 9):
        plot3d(facets, data, values_dim)
    elif (facets.shape[1] == 4):
        plot2d(facets, data, values_dim)
    else:
        raise Exception("Facets are not 2D or 3D! Corrupt file.")


if __name__ == "__main__":
    advice = "Usage is: python py/data_plotter.py filename column" +\
             "\n where the column specifies the column of the values" +\
             " dataspace to display"
    if len(sys.argv) < 3:
        print(advice)
        sys.exit()
    filename = sys.argv[1]
    try:
        values_dim = int(sys.argv[2])
    except:
        print(advice)
    main(sys.argv[1], int(sys.argv[2]))
