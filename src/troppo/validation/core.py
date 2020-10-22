from cobamp.wrappers.core import ConstraintBasedModelSimulator

from numpy import nan


def default_post_processing_func(x):
    status, obj, values = x
    if not status:
        obj = nan
        values = {k:nan for k in values.keys()}
    return values, obj


class ContextSpecificModelSimulator(object):
    def __init__(self, model_simulator: ConstraintBasedModelSimulator, flux_data: dict, scenarios: dict,
                 post_process=None):
        self.model_simulator = model_simulator
        self.data = flux_data
        self.scenarios = scenarios
        self.post_process = default_post_processing_func if post_process is not None else post_process

    def __simulate(self, contexts: dict, simulation_function, objective_coefficients: dict, minimize: bool, **kwargs):
        context_names, context_bounds = zip(contexts.items())

        full_cnames, full_cbounds = [], []
        for sc_name, sc_bounds in self.scenarios:
            full_cnames.extend([(sc_name, c_name) for c_name, c_bounds in contexts.items()])
            context_bounds_scen = [{k:v for k,v in d.items()} for d in context_bounds]
            for d in context_bounds_scen: d.update(sc_bounds)
            full_cbounds.extend(context_bounds_scen)

        simulations = zip(*self.model_simulator.batch_simulate(func=simulation_function,
                                            bound_changes=full_cbounds,
                                            objective_coefficients=[objective_coefficients],
                                            minimize=[minimize],
                                            **kwargs))

        return dict(zip(full_cnames,[self.post_process(k) for k in simulations]))