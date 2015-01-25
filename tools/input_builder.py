import numpy as np
from mako.template import Template

class Element(object):
    def __init__(self, pts, bc_type, bc, refine):
        self.pts = np.array(pts).astype(np.float32).tolist()
        self.bc_type = str(bc_type)
        self.bc = np.array(bc).astype(np.float32).tolist()
        self.refine = int(refine)

def displacement_edge(end_pts, refine, fnc):
    n = 2 ** refine
    x_vals = np.linspace(end_pts[0][0], end_pts[1][0], n)
    y_vals = np.linspace(end_pts[0][1], end_pts[1][1], n)
    ux, uy = fnc(x_vals, y_vals)

    es = []
    for i in range(n - 1):
        es.append(Element(
            [[x_vals[i], y_vals[i]], [x_vals[i + 1], y_vals[i + 1]]],
            "displacement",
            [[ux[i], uy[i]], [ux[i + 1], uy[i + 1]]],
            0
        ))

    return es

def exec_template(filename, **params):
    file_template = """
    {
        "shear_modulus": ${G},
        "poisson_ratio": ${mu},
        "elements": [
        % for e in es:
            {
                "pts": ${e.pts},
                "bc_type": "${e.bc_type}",
                "bc": ${e.bc},
                "refine": ${e.refine}
            }
            % if loop.index != len(es) - 1:
            ,
            % endif
        % endfor
        ]
    }
    """

    text = Template(file_template).render(**params)
    with open(filename, 'w') as file:
        file.write(text)
