import bie_spec
from log_tools import log_exceptions
from interface import Executor
import constraints
import system
import meshing

import copy
import numpy as np
import logging
logger = logging.getLogger(__name__)

@log_exceptions
def adaptive_execute(dim, elements, input_params, callback = lambda x, y: None):
    log_begin_adaptive(len(elements), input_params)
    executor = AdaptiveExecutor(dim, elements, input_params, callback)
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
    def __init__(self, dim, elements, input_params, callback):
        self.arguments = (dim, elements, input_params)
        self.error_threshold = input_params['error_threshold']
        self.max_iters = input_params['max_iters']
        self.refine_fraction = input_params['refine_fraction']
        self.callback = callback

    def run(self):
        problem = copy.copy(self.arguments)
        for i in range(self.max_iters):
            executor = Executor(*problem)
            print(sum([v.n_facets() for v in executor.meshes.values()]))
            result = executor.run()
            residual = self.calculate_residual(executor, result.soln)
            self.callback(executor, result)
            local_error, global_error = self.error_estimate(executor, residual)
            print(global_error)
            if global_error < self.error_threshold:
                break
            which_to_refine = self.choose_refines(executor, local_error)
            new_es = self.refine(executor.element_lists, which_to_refine)
            problem = (problem[0], new_es, problem[2])
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
            # facet_error *= facet_lengths
            local_error[m] = facet_error
        all_errors = np.concatenate(local_error.values())
        global_error = np.sqrt(np.mean(all_errors ** 2))
        return local_error, global_error

    def choose_refines(self, executor, local_error):
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
        return separated_by_mesh

    def refine(self, element_lists, which_to_refine):
        es = []
        for m in bie_spec.mesh_types:
            es.extend(self.refine_mesh(
                element_lists[m], {k:True for k in which_to_refine[m]}
            ))
        return es

    def refine_mesh(self, es, which_to_refine):
        out = []
        for i in range(len(es)):
            if which_to_refine.get(i, False) is False:
                out.append(es[i])
                continue
            midpt = (es[i]['pts'][0] + es[i]['pts'][1]) / 2.0
            out.append(dict(
                type = es[i]['type'],
                pts = [es[i]['pts'][0], midpt],
                constraints = [
                   constraints.refine_constraint(c, 0) for c in es[i]['constraints']
                ]
            ))
            out.append(dict(
                type = es[i]['type'],
                pts = [midpt, es[i]['pts'][1]],
                constraints = [
                   constraints.refine_constraint(c, 1) for c in es[i]['constraints']
                ]
            ))
        return out

    def log_adaptive_step(self, step_idx, error):
        logger.info(
            'Adaptive refinement step #' + str(step_idx) + ' complete ' +
            'with error: ' + str(error)
        )
