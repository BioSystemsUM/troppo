from cobamp.wrappers.core import ConstraintBasedModelSimulator
from cobamp.utilities.hash import get_unique_dicts

from numpy import nan


def default_post_processing_func(x):
    status, obj, values = x
    if not status:
        obj = nan
        values = {k:nan for k in values.keys()}
    return values, obj


class ContextSpecificModelSimulator(object):
    def __init__(self, model_simulator: ConstraintBasedModelSimulator, scenarios=None,
                 post_process=None):
        self.model_simulator = model_simulator
        self.scenarios = scenarios if scenarios is not None else {'default':{}}
        self.post_process = default_post_processing_func if post_process is not None else post_process


    def simulate(self, contexts: dict, simulation_function, objective_coefficients: dict, minimize: bool, **kwargs):
        context_names, context_bounds = zip(contexts.items())

        full_cnames, full_cbounds = [], []
        for sc_name, sc_bounds in self.scenarios:

            full_cnames.extend([(sc_name, c_name) for c_name, c_bounds in contexts.items()])
            context_bounds_scen = [{k:v for k,v in d.items()} for d in context_bounds]
            for d in context_bounds_scen: d.update(sc_bounds)
            full_cbounds.extend(context_bounds_scen)

        dict_cache, dict_indices = get_unique_dicts(full_cbounds)
        dnames, dobjs = zip(dict_cache.items())

        simulations = zip(*self.model_simulator.batch_simulate(func=simulation_function,
                                            bound_changes=dobjs,
                                            objective_coefficients=[objective_coefficients],
                                            minimize=[minimize],
                                            **kwargs))

        pre_result = dict(zip(dnames,[self.post_process(k) for k in simulations]))

        result_list = [pre_result[k] for k in dict_indices]
        return dict(zip(full_cnames, result_list))
