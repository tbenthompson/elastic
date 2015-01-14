import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import sys

def main(filename):
    f = h5py.File(filename)
    locs = f['locations']
    data = f['values']

    skip = 1
    x = locs[:, 0]
    y = locs[:, 1]
    ux = data[:, 0]
    uy = data[:, 1]
    plt.quiver(x, y, ux, uy, pivot = 'tail')
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.show()

if __name__ == "__main__":
    advice = "Usage is: python py/data_plotter.py filename"
    if len(sys.argv) < 2:
        print(advice)
        sys.exit()
    filename = sys.argv[1]
    main(filename)
