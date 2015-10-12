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

def evaluate_interior(dispatcher, soln, pts, normals, bie):
    dim = pts.shape[1]
    result = np.zeros((dim, pts.shape[0]))

    #TODO: This chunking is only necessary because dense matrices are being used
    # move to using FMM.
    chunk_size = 1000
    chunks = int(np.ceil(pts.shape[0] / float(chunk_size)))
    for i in range(chunks):
        start_idx = i * chunk_size
        past_end_idx = start_idx + min(pts.shape[0] - start_idx + 1, chunk_size)
        chunk_pts = pts[start_idx:past_end_idx]
        chunk_normals = normals[start_idx:past_end_idx]
        for t in bie['terms']:
            op = dispatcher.compute_interior(t, chunk_pts, chunk_normals)
            f = op.select_input_field(soln)
            if f is None:
                return None
            # Negate to shift to the RHS.
            # The integral equation is set up like u(x) + Integrals = 0.
            # We want u(x) = -Integrals
            chunk_result = -op.apply(f)
            for d in range(dim):
                dim_start_idx = d * chunk_pts.shape[0]
                dim_end_idx = (d + 1) * chunk_pts.shape[0]
                result[d, start_idx:past_end_idx] +=\
                    chunk_result[dim_start_idx:dim_end_idx]

    result *= bie['mass_term']['multiplier']

    #TODO: Is this necessary? Just leave it as a big numpy array.
    out = []
    for d in range(dim):
        out.append(result[d, :])
    return out


def scale(unknowns, scalings, params, inverse):
    for u, values in unknowns.iteritems():
        factor = scalings[u](params)
        if inverse:
            factor = 1.0 / factor
        for d in range(len(values)):
            values[d] *= factor

def add_constant_fields(unknowns, zeros):
    scale_factor = 1.0
    if zeros:
        scale_factor = 0.0
    dim = len(unknowns[('continuous', 'displacement')])
    n_dofs_continuous = len(unknowns[('continuous', 'displacement')][0])
    unknowns[('continuous', 'ones')] = [
        scale_factor * np.ones(n_dofs_continuous) for d in range(dim)
    ]
