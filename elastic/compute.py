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
        if self.spec['function'] == 'ones':
            f = [np.ones(self.internal.n_cols() / dim) for d in range(dim)]
        return f

    def data(self):
        if type(self.internal) is tbempy._tbempy.DenseOperator:
            return self.internal.data()
        elif type(self.internal) is tbempy._tbempy.SparseOperator:
            return self.internal.to_dense().data()
        return None

    def to_numpy_matrix(self):
        rows = self.internal.n_rows()
        cols = self.internal.n_cols()
        return self.spec['multiplier'] * self.data().reshape((rows, cols))

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
        result = self.evaluator.boundary(obs_mesh, src_mesh, kernel)
        return Op(result, op_spec, False)

    def compute_interior(self, op_spec, pts, normals):
        src_mesh = self.meshes[op_spec['src_mesh']]
        kernel = self.kernels[op_spec['kernel']]
        result = self.evaluator.interior(pts, normals, src_mesh, kernel)
        return Op(result, op_spec, False)

    def compute_mass(self, mass_spec):
        obs_mesh = self.meshes[mass_spec['src_mesh']]
        return Op(self.evaluator.mass(obs_mesh), mass_spec, True)

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
