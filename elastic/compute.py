import numpy as np
import tbempy

class Op(object):
    def __init__(self, op, spec, is_mass):
        self.internal = op
        self.spec = spec
        self.is_mass = is_mass

    def apply(self, field):
        computed = self.internal.apply(np.concatenate(field))
        return computed * self.spec['multiplier']

    def input_type(self):
        return (self.spec['src_mesh'], self.spec['function'])

    def select_input_field(self, fields):
        dim = len(fields.values()[0])
        f = fields.get(self.input_type(), None)
        return f

    def to_dense(self):
        if type(self.internal) is tbempy._tbempy.DenseOperator:
            return self.internal
        elif type(self.internal) is tbempy._tbempy.SparseOperator:
            return self.internal.to_dense()

    def get_multiplier(self):
        return self.spec['multiplier']

class IntegralEvaluator(object):
    def mass(self, obs_mesh):
        return self.tbem.mass_operator_tensor(obs_mesh, self.params['obs_far_order'])

    def make_integrator(self, kernel):
        return self.tbem.make_sinh_integrator(self.params['sinh_order'],
            self.params['obs_near_order'], self.params['obs_far_order'],
            self.params['src_far_order'], self.params['singular_steps'],
            self.params['far_threshold'], kernel
        )

class FMMIntegralEvaluator(IntegralEvaluator):
    def __init__(self, tbem, params, all_mesh):
        self.tbem = tbem
        self.params = params
        self.fmm_config = tbem.FMMConfig(0.3, params['fmm_order'], 250, 0.1, True)
        self.all_mesh = all_mesh

    def boundary(self, obs_mesh, src_mesh, kernel):
        mthd = self.make_integrator(kernel)
        return self.tbem.boundary_operator(
            obs_mesh, src_mesh, mthd, self.fmm_config, self.all_mesh
        )

    def interior(self, pts, normals, src_mesh, kernel):
        mthd = self.make_integrator(kernel)
        #TODO: Switch to non-dense once available!
        return self.tbem.dense_interior_operator(
            pts, normals, src_mesh, mthd, self.all_mesh
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
    def __init__(self, mesh_provider, kernels, evaluator):
        self.mesh_provider = mesh_provider
        self.kernels = kernels
        self.evaluator = evaluator

    def compute_boundary(self, op_spec):
        obs_mesh = self.mesh_provider.get_obs_mesh(op_spec)
        src_mesh = self.mesh_provider.get_src_mesh(op_spec)
        kernel = self.kernels[op_spec['kernel']]
        result = self.evaluator.boundary(obs_mesh, src_mesh, kernel)
        final_op = self.mesh_provider.distribute_zeros(op_spec, result)
        out = Op(final_op, op_spec, False)
        return out

    def compute_interior(self, op_spec, pts, normals):
        src_mesh = self.mesh_provider[op_spec['src_mesh']]
        kernel = self.kernels[op_spec['kernel']]
        result = self.evaluator.interior(pts, normals, src_mesh, kernel)
        return Op(result, op_spec, False)

    def compute_mass(self, mass_spec):
        obs_mesh = self.mesh_provider.get_src_mesh(mass_spec)
        return Op(self.evaluator.mass(obs_mesh), mass_spec, True)

class BIE(object):
    def __init__(self, terms, spec, unknowns_to_knowns):
        self.terms = terms
        self.spec = spec
        self.unknowns_to_knowns = unknowns_to_knowns

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
        return (self.spec['obs_mesh'], self.spec['unknown_field'])
