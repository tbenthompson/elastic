from compute import BIE
import bie_spec
import numpy as np

def form_linear_systems(bies, dispatcher):
    systems = []
    for spec in bies:
        integrals = []
        for term in spec['terms']:
            integrals.append(dispatcher.compute_boundary(term))
        integrals.append(dispatcher.compute_mass(spec['mass_term']))
        systems.append(BIE(integrals, spec, bie_spec.unknowns_to_knowns))
    return systems

def evaluate_linear_systems(bies, fields):
    dim = len(fields.values()[0])
    result = dict()
    for b in bies:
        assert(b.output_type() not in result)
        rows = b.evaluate(fields)
        assert(rows is not None)
        result[b.output_type()] = split_into_components(dim, rows)
    return result

def split_into_components(dim, field):
    result = []
    n_dofs_per_components = field.shape[0] / dim
    for d in range(dim):
        start_dof = n_dofs_per_components * d
        end_dof = n_dofs_per_components * (d + 1)
        result.append(field[start_dof:end_dof])
    return result

def evaluate_interior(dispatcher, soln, pts, normals, terms):
    dim = pts.shape[1]
    result = np.zeros(pts.shape[0] * dim)
    for t in terms:
        op = dispatcher.compute_interior(t, pts, normals)
        f = op.select_input_field(soln)
        if f is None:
            return None
        # Negate to shift to the RHS.
        # The integral equation is set up like u(x) + Integrals = 0.
        # We want u(x) = -Integrals
        result -= op.apply(f)
    out = []
    for d in range(dim):
        start_idx = d * pts.shape[0]
        end_idx = (d + 1) * pts.shape[0]
        out.append(result[start_idx:end_idx])
    return out

def scale(unknowns, scalings, params, inverse):
    for u, values in unknowns.iteritems():
        factor = scalings[u](params)
        if inverse:
            factor = 1.0 / factor
        for d in range(len(values)):
            values[d] *= factor
