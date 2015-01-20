import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import sys

def main(filenamex, filenamey):
    fx = h5py.File(filenamex)
    fy = h5py.File(filenamey)
    locs = fx['locations']
    datax = fx['values'][:,0]
    datay = fy['values'][:,0]

    skip = 1
    x = locs[:, 0]
    y = locs[:, 1]
    plt.quiver(x, y, datax, datay, pivot = 'tail')
    plt.plot([-0.5, 0.5], [-0.5, 0.5], 'b')
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.show()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
