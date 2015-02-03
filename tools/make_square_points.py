from input_builder import exec_template
import sys
import numpy as np

def main():
    filename = sys.argv[1]
    x_guide = [float(v) for v in sys.argv[2].split(',')]
    y_guide = [float(v) for v in sys.argv[3].split(',')]

    x_vals = np.linspace(*x_guide)
    y_vals = np.linspace(*y_guide)
    x_pts, y_pts = np.meshgrid(x_vals, y_vals)
    ps = zip(x_pts.flatten(), y_pts.flatten())

    file_template = """
    [
        % for p in ps:
        [${p[0]}, ${p[1]}]
        % if loop.index != len(ps) - 1:
        ,
        % endif
        % endfor
    ]
    """
    exec_template(file_template, filename, ps = ps)


if __name__ == "__main__":
    main()
