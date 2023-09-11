import numpy as np
from collections import OrderedDict
from numpy import ndarray, array

from cobamp.core.models import ConstraintBasedModel
from cobamp.core.optimization import Solution
from troppo.methods.base import ContextSpecificModelReconstructionAlgorithm, PropertiesReconstruction


class GIMMEModel(ConstraintBasedModel):
    """
    A GIMME model is a model that is used to reconstruct a context-specific model using the GIMME algorithm.

    Parameters:
    -----------
    cbmodel: ConstraintBasedModel
        The template model that is used to reconstruct the context-specific model.
    solver: str or None
        The solver that is used to solve the optimization problem.

    Attributes:
    -----------
    cbmodel: ConstraintBasedModel
        The template model that is used to reconstruct the context-specific model.
    mapping: dict
        A dictionary that maps the reactions of the template model to the reactions of the GIMME model.

    """
    def __init__(self, cbmodel: ConstraintBasedModel, solver: str or None = None):
        self.cbmodel = cbmodel
        if not self.cbmodel.model:
            self.cbmodel.initialize_optimizer()

        irrev_model, self.mapping = cbmodel.make_irreversible()

        S = irrev_model.get_stoichiometric_matrix()
        bounds = irrev_model.bounds
        super().__init__(S, bounds, irrev_model.reaction_names, irrev_model.metabolite_names, solver=solver,
                         optimizer=True)

    def __adjust_objective_to_irreversible(self, objective_dict):
        obj_dict = {}
        for k, v in objective_dict.items():
            irrev_map = self.mapping[self.cbmodel.decode_index(k, 'reaction')]
            if isinstance(irrev_map, (list, tuple)):
                for i in irrev_map:
                    obj_dict[i] = v
            else:
                obj_dict[irrev_map] = v
        return obj_dict

    def __adjust_expression_vector_to_irreversible(self, exp_vector):
        exp_vector_n = np.zeros(len(self.reaction_names), )
        for rxn, val in enumerate(exp_vector):
            rxmap = self.mapping[rxn]
            if isinstance(rxmap, tuple):
                exp_vector_n[rxmap[0]] = exp_vector_n[rxmap[1]] = val
            else:
                exp_vector_n[rxmap] = val
        return exp_vector_n

    def optimize_gimme(self, exp_vector: list, objectives: list or tuple, obj_frac: list or tuple or float = 0.9,
                       flux_thres: float = None):
        """
        Optimize the GIMME model.

        Parameters
        ----------
        exp_vector: list
            A list of expression values for each reaction in the GIMME model.
        objectives: list or tuple
            A list of dictionaries that define the objectives of the GIMME model.
        obj_frac: list or tuple or float
            A list of fractions that define the lower bounds of the objectives. If a float is given, the same fraction
            is used for all objectives.
        flux_thres: float
            A threshold that defines the minimum flux that is allowed for each reaction in the GIMME model.

        Returns
        -------
        GIMMESolution:
            The solution of the GIMME model.
        """
        N = len(self.cbmodel.reaction_names)
        objectives_irr = [self.__adjust_objective_to_irreversible(obj) for obj in objectives]
        exp_vector_irr = self.__adjust_expression_vector_to_irreversible(exp_vector)

        def find_objective_value(obj):
            self.cbmodel.set_objective(obj, False)
            return self.cbmodel.optimize().objective_value()

        objective_values = list(map(find_objective_value, objectives_irr))

        gimme_model_objective = array(
            [flux_thres - exp_vector_irr[i] if -1 < exp_vector_irr[i] < flux_thres else 0 for i in range(N)])

        objective_lbs = np.zeros(len(self.reaction_names))
        for ov, obj in zip(objective_values, objectives_irr):
            for rx, v in obj.items():
                objective_lbs[rx] = v * ov * obj_frac

        objective_ids = np.nonzero(objective_lbs)[0]
        lbs_id = objective_lbs[objective_ids]
        for idx, lb in zip(objective_ids, lbs_id):
            self.set_reaction_bounds(idx, lb=lb, temporary=True)

        self.set_objective(gimme_model_objective, True)
        sol = self.optimize()
        self.revert_to_original_bounds()
        return GIMMESolution(sol, exp_vector, self.cbmodel.reaction_names, self.mapping)


class GIMMESolution(Solution):
    """
    A GIMME solution is a solution of a GIMME model.

    Parameters:
    -----------
    sol: Solution
        The solution of the GIMME model.
    exp_vector: list or ndarray
        Expression vector
    var_names: list
        List of variable names
    mapping: dict
        A dictionary that maps the reactions of the template model to the reactions of the GIMME model.

    """
    def __init__(self, sol, exp_vector, var_names, mapping=None):
        self.exp_vector = exp_vector
        gimme_solution = sol.x()
        if mapping:
            gimme_solution = [max(gimme_solution[array(new)]) if isinstance(new, (tuple, list)) else gimme_solution[new]
                              for orig, new
                              in mapping.items()]
        super().__init__(
            value_map=OrderedDict([(k, v) for k, v in zip(var_names, gimme_solution)]),
            status=sol.status(),
            objective_value=sol.objective_value()
        )

    def get_reaction_activity(self, flux_threshold: float):
        """
        Get the reaction activity of the GIMME solution.

        Parameters
        ----------
        flux_threshold: float
            The flux threshold that is used to determine the reaction activity.

        Returns
        -------
        reaction_index: list
            List of reaction indices from active reactions.
        """
        gimme_fluxes = array([kv[1] for i, kv in enumerate(self.var_values().items())])
        activity = np.zeros(gimme_fluxes.shape)
        ones = (np.array(self.exp_vector) > flux_threshold) | (self.exp_vector == -1)
        twos = gimme_fluxes > 0
        activity[ones] = 1
        activity[twos & ~ones] = 2

        reaction_index = [idx for idx, val in enumerate(activity) if val != 0]
        return reaction_index


class GIMMEProperties(PropertiesReconstruction):
    """
    Properties for GIMME

    Parameters:
    -----------
    exp_vector: list or ndarray
        Expression vector
    objectives: list or tuple
        List of objective vectors
    obj_frac: list or tuple or float
        Fraction of the objective vector to be used
    preprocess: bool
        Preprocess the model
    flux_threshold: float
        Flux threshold
    solver: str
        Solver to be used
    reaction_ids: list
        List of reaction ids
    metabolite_ids: list
        List of metabolite ids

    """
    def __init__(self, exp_vector: list, objectives: list or tuple, obj_frac: list or tuple or float = 0.9,
                 preprocess: bool = False, flux_threshold: float = None, solver: str = None, reaction_ids: list = None,
                 metabolite_ids: list = None):
        new_mandatory = {
            'exp_vector': lambda x: isinstance(x, list) and len(x) > 0 or isinstance(x, ndarray),
            'preprocess': lambda x: isinstance(x, bool) or x is None,
            'objectives': lambda x: type(x) in [list, tuple, ndarray],
            'reaction_ids': lambda x: isinstance(x, list) and len(x) > 0 or isinstance(x, ndarray),
            'metabolite_ids': lambda x: isinstance(x, list) and len(x) > 0 or isinstance(x, ndarray)}

        new_optional = {'obj_frac': lambda x: type(x) in [ndarray, list, tuple, float],
                        'flux_threshold': lambda x: isinstance(x, float) or x is None,
                        'solver': lambda x: isinstance(x, str) or x is None}
        super().__init__()

        self.add_new_properties(new_mandatory, new_optional)

        self['objectives'] = objectives
        self['exp_vector'] = exp_vector
        self['solver'] = solver
        self['reaction_ids'] = reaction_ids
        self['metabolite_ids'] = metabolite_ids
        self['obj_frac'] = obj_frac if isinstance(obj_frac, ndarray) else array([obj_frac] * len(objectives))
        self['preprocess'] = True if preprocess else False
        self['flux_threshold'] = 1e-4 if flux_threshold is None else flux_threshold

    @staticmethod
    def from_integrated_scores(scores: list, **kwargs):
        """
        Create GIMMEProperties from integrated scores

        Parameters
        ----------
        scores: list or ndarray
            Integrated scores
        kwargs: dict
            Additional arguments

        Returns
        -------
        GIMMEProperties

        """
        return GIMMEProperties(exp_vector=scores, **{k: v for k, v in kwargs.items() if 'exp_vector' not in k})


class GIMME(ContextSpecificModelReconstructionAlgorithm):
    """
    GIMME algorithm

    Parameters
    ----------
    S: list or ndarray
        Stoichiometric matrix
    lb: list or ndarray
        Lower bounds
    ub: list or ndarray
        Upper bounds
    properties: GIMMEProperties
        GIMME properties

    Attributes
    ----------
    S: ndarray
        Stoichiometric matrix
    lb: ndarray
        Lower bounds
    ub: ndarray
        Upper bounds
    properties: GIMMEProperties
        GIMME properties
    sol: GIMMESolution
        GIMME solution
    gm: GIMMEModel
        GIMME model
    """
    properties_class = GIMMEProperties

    def __init__(self, S: list, lb: list, ub: list, properties: GIMMEProperties):
        super().__init__(S, lb, ub, properties)
        self.S = np.array(S)
        self.lb, self.ub = np.array(lb), np.array(ub)
        self.properties = properties
        self.model = GIMMEModel
        self.sol = None
        cbm = ConstraintBasedModel(S, list(zip(lb, ub)), reaction_names=self.properties['reaction_ids'],
                                   metabolite_names=self.properties['metabolite_ids'])
        self.gm = GIMMEModel(cbm, self.properties['solver'])

    def run(self):
        """
        Run GIMME algorithm

        Returns
        -------
        list: List with index of active reactions.

        """
        sol = self.gm.optimize_gimme(
            exp_vector=self.properties['exp_vector'],
            objectives=self.properties['objectives'],
            obj_frac=self.properties['obj_frac'],
            flux_thres=self.properties['flux_threshold']
        )
        self.sol = sol
        return sol.get_reaction_activity(self.properties['flux_threshold'])
