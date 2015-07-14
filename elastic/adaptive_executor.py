import bie_spec
from exceptions import log_exceptions
from interface import Executor, create_input
import system
import meshing

import numpy as np
import logging
logger = logging.getLogger(__name__)

@log_exceptions
def adaptive_execute(dim, elements, input_params):
    log_begin_adaptive(len(elements), input_params)
    executor = AdaptiveExecutor(dim, elements, input_params)
    return executor.run()

def log_begin_adaptive(n_elements, input_params):
    logger.info('')
    logger.info(
        'Beginning adaptive elastic calculation with ' +
        str(n_elements) +
        'starting elements and with with params: ' +
        str(input_params)
    )

class AdaptiveExecutor(object):
    def __init__(self, dim, elements, input_params):
        self.arguments = (dim, elements, input_params)
        self.initial_problem = create_input(*self.arguments)
        self.error_threshold = input_params['error_threshold']
        self.max_iters = input_params['max_iters']
        self.refine_fraction = input_params['refine_fraction']

    def run(self):
        problem = self.initial_problem
        for i in range(self.max_iters):
            print(sum([v.n_facets() for v in problem[2].values()]))
            executor = Executor(*problem)
            result = executor.run()
            residual = self.calculate_residual(executor, result.soln)
            local_error, global_error = self.error_estimate(executor, residual)
            print(global_error)
            if global_error < self.error_threshold:
                break
            new_meshes = self.refine(executor, local_error)
            meshing.postprocess_meshes(executor.tbem, new_meshes)
            problem = (problem[0], problem[1], new_meshes)
            self.log_adaptive_step(i, global_error)
        return result

    def calculate_residual(self, executor, solution):
        bies = bie_spec.get_residual_BIEs(executor.params)
        residual = system.evaluate_linear_systems(executor.assemble(bies), solution)
        system.scale(residual, bie_spec.integral_scaling, executor.params, True)
        return residual

    def error_estimate(self, executor, residual):
        local_error = dict()
        for m in bie_spec.mesh_types:
            fs = executor.meshes[m].facets
            facet_lengths = np.linalg.norm(fs[:, 0, :] - fs[:, 1, :], axis = 1)
            facet_error = np.zeros_like(facet_lengths)
            for mesh_name, field_name in residual:
                if mesh_name != m:
                    continue
                for d in range(2):
                    dof_error = np.abs(residual[(mesh_name, field_name)][d])
                    facet_error += (dof_error[::2] + dof_error[1::2])
            facet_error *= facet_lengths
            local_error[m] = facet_error
        all_errors = np.concatenate(local_error.values())
        global_error = np.sqrt(np.mean(all_errors ** 2))
        return local_error, global_error

    def refine(self, executor, local_error):
        tracers = []
        for m in bie_spec.mesh_types:
            fs = executor.meshes[m]
            tracers.extend([(m, i) for i in range(fs.n_facets())])
        worst_to_best = sorted(tracers, key = lambda k: local_error[k[0]][k[1]])
        n_facets = len(worst_to_best)
        n_refine = int(np.ceil(n_facets * self.refine_fraction))
        refine_me = worst_to_best[-n_refine:]
        #TODO: log number of facets being refined
        separated_by_mesh = {m:[] for m in bie_spec.mesh_types}
        for item in refine_me:
            separated_by_mesh[item[0]].append(item[1])
        out = dict()
        for m in bie_spec.mesh_types:
            out[m] = executor.meshes[m].refine(separated_by_mesh[m])
        return out

    def log_adaptive_step(self, step_idx, error):
        logger.info(
            'Adaptive refinement step #' + str(step_idx) + ' complete ' +
            'with error: ' + str(error)
        )
