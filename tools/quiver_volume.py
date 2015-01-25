import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import sys

def main(filename):
    f = h5py.File(filename)
    locs = f['locations']
    datax = f['values0'][:,0]
    datay = f['values1'][:,0]

    skip = 1
    x = locs[:, 0]
    y = locs[:, 1]
    plt.quiver(x, y, datax, datay, pivot = 'tail')
    # plt.plot([-0.5, 0.5], [-0.5, 0.5], 'b')
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.show()

if __name__ == "__main__":
    main(sys.argv[1])
