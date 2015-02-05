from input_builder import points_grid
import sys
import numpy as np

def main():
    filename = sys.argv[1]
    x_range = [float(v) for v in sys.argv[2].split(',')]
    y_range = [float(v) for v in sys.argv[3].split(',')]
    points_grid(x_range, y_range, filename)

if __name__ == "__main__":
    main()
