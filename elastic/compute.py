import numpy as np

class Op(object):
    def __init__(self, op, spec):
        self.internal = op
        self.spec = spec

    def function(self):
        return self.spec['function']

    def src_mesh(self):
        return self.spec['src_mesh']

    def apply(self, field):
        computed = self.internal.apply(np.concatenate(field))
        return -computed * self.spec['multiplier']

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
    def __init__(self, terms):
        self.terms = terms

    def evaluate(self, fields):
        n_rows = self.terms[0].internal.n_rows()
        result = np.zeros(n_rows)
        for t in self.terms:
            f = fields.get((t.src_mesh(), t.function()), None)
            if f is None:
                print(t.spec)
                return None
            result += t.apply(f)
        return result

#TODO: These don't make sense here.
#To an "evaluate" module? Same place as scale columns/rows?
def form_linear_systems(bies, dispatcher):
    systems = []
    for b in bies:
        integrals = []
        for term in b['terms']:
            integrals.append(dispatcher.compute_boundary(term))
        integrals.append(dispatcher.compute_mass(b['mass_term']))
        systems.append(BIE(integrals))
    return systems

def evaluate_linear_systems(bies, fields):
    result = []
    for b in bies:
        evaluated = b.evaluate(fields)
        assert(evaluated is not None)
        result.append(evaluated)
    return result

# def interior_eval(tbem, meshes, kernels, params, soln, pts, normals, terms):
#     for k, v in soln.iteritems():
#         fields[k] = v
#     result = np.zeros(pts.shape[0] * tbem.dim)
#     for t in terms:
#         op = compute_interior_integral(
#             tbem, meshes, kernels, params, t, pts, normals, fields
#         )
#         result += apply_op(op, fields)
#     out = []
#     for d in range(tbem.dim):
#         start_idx = d * pts.shape[0]
#         end_idx = (d + 1) * pts.shape[0]
#         out.append(result[start_idx:end_idx])
#     return out
