import numpy as np

class Op(object):
    def __init__(self, op, spec):
        self.internal = op
        self.spec = spec

    def apply(self, field):
        computed = self.internal.apply(np.concatenate(field))
        return computed * self.spec['multiplier']

    def select_input_field(self, fields):
        dim = len(fields.values()[0])
        f = fields.get((self.spec['src_mesh'], self.spec['function']), None)
        if self.spec['function'] == 'ones':
            f = [np.ones(self.internal.n_cols() / dim) for d in range(dim)]
        return f

class IntegralEvaluator(object):
    def mass(self, obs_mesh):
        return self.tbem.mass_operator_tensor(obs_mesh, self.params['obs_order'])

    def make_integrator(self, kernel):
        return self.tbem.make_sinh_integrator(
            self.params['sinh_order'], self.params['obs_order'],
            self.params['singular_steps'], self.params['far_threshold'],
            kernel
        )

class FMMIntegralEvaluator(IntegralEvaluator):
    def __init__(self, tbem, params, all_mesh):
        self.tbem = tbem
        self.params = params
        self.fmm_config = tbem.FMMConfig(0.3, params['fmm_order'], 250, 0.1, True)
        self.all_mesh = all_mesh

    def boundary(self, obs_mesh, src_mesh, kernel):
        mthd = self.make_integrator(kernel)
        return tbem.boundary_operator(
            obs_mesh, src_mesh, mthd, self.fmm_config, self.all_mesh
        )

    def interior(self, pts, normals, src_mesh, kernel):
        mthd = self.make_integrator(kernel)
        return tbem.interior_operator(
            pts, normals, src_mesh, mthd, self.fmm_config, self.all_mesh
        )

class DenseIntegralEvaluator(IntegralEvaluator):
    def __init__(self, tbem, params, all_mesh):
        self.tbem = tbem
        self.params = params
        self.all_mesh = all_mesh

    def boundary(self, obs_mesh, src_mesh, kernel):
        mthd = self.make_integrator(kernel)
        return self.tbem.dense_boundary_operator(
            obs_mesh, src_mesh, mthd, self.all_mesh
        )

    def interior(self, pts, normals, src_mesh, kernel):
        mthd = self.make_integrator(kernel)
        return self.tbem.dense_interior_operator(
            pts, normals, src_mesh, mthd, self.all_mesh
        )

class IntegralDispatcher(object):
    def __init__(self, meshes, kernels, evaluator):
        self.meshes = meshes
        self.kernels = kernels
        self.evaluator = evaluator

    def compute_boundary(self, op_spec):
        obs_mesh = self.meshes[op_spec['obs_mesh']]
        src_mesh = self.meshes[op_spec['src_mesh']]
        kernel = self.kernels[op_spec['kernel']]
        return Op(self.evaluator.boundary(obs_mesh, src_mesh, kernel), op_spec)

    def compute_interior(self, op_spec, pts, normals):
        src_mesh = self.meshes[op_spec['src_mesh']]
        kernel = self.kernels[op_spec['kernel']]
        return Op(self.evaluator.interior(pts, normals, src_mesh, kernel), op_spec)

    def compute_mass(self, mass_spec):
        obs_mesh = self.meshes[mass_spec['src_mesh']]
        return Op(self.evaluator.mass(obs_mesh), mass_spec)

class BIE(object):
    def __init__(self, terms, spec):
        self.terms = terms
        self.spec = spec

    def evaluate(self, fields):
        dim = len(fields.values()[0])
        n_rows = self.terms[0].internal.n_rows()
        result = np.zeros(n_rows)
        for t in self.terms:
            f = t.select_input_field(fields)
            if f is None:
                return None
            result += t.apply(f)
        return result

    def output_type(self):
        return (self.spec['obs_mesh'], self.spec['mass_term']['function'])

#TODO: These don't make sense here.
#To an "evaluate" module? Same place as scale columns/rows?
#Or a "system" module?
def form_linear_systems(bies, dispatcher):
    systems = []
    for spec in bies:
        integrals = []
        for term in spec['terms']:
            integrals.append(dispatcher.compute_boundary(term))
        integrals.append(dispatcher.compute_mass(spec['mass_term']))
        systems.append(BIE(integrals, spec))
    return systems

def split_into_components(dim, field):
    result = []
    n_dofs_per_components = field.shape[0] / dim
    for d in range(dim):
        start_dof = n_dofs_per_components * d
        end_dof = n_dofs_per_components * (d + 1)
        result.append(field[start_dof:end_dof])
    return result

def evaluate_linear_systems(bies, fields):
    dim = len(fields.values()[0])
    result = dict()
    for b in bies:
        evaluated = b.evaluate(fields)
        assert(evaluated is not None)
        result[b.output_type()] = split_into_components(dim, evaluated)
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
