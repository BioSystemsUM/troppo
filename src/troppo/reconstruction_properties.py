from cobamp.utilities.property_management import PropertyDictionary
from cobamp.wrappers.external_wrappers import model_readers

from numbers import Number
from troppo.methods_reconstruction import MethodsReconstruction
from numpy import ndarray, array
from multiprocessing import cpu_count


# Help functions

def is_list(x):
	return type(x) in [list, tuple] and len(x) > 0 or isinstance(x, ndarray) and x.size > 0


def is_list_else_empty(x):
	return type(x) in [list, tuple] or isinstance(x, ndarray)


def if_none_return_list(x):
	if x is None:
		return array([])
	else:
		return x


class PropertiesReconstruction(PropertyDictionary):
	def __init__(self):
		self.base_mandatory = {'solver': str, 'template_model': model_readers, 'method': MethodsReconstruction,
							   'omics_type': 'omics'}
		self.base_optional = {'env_conditions': dict, 'flux_threshold': float}
		super().__init__(self.base_mandatory, self.base_optional)


class FastcoreProperties(PropertiesReconstruction):
	def __init__(self, core, flux_threshold=1e-4, solver=None):
		new_mandatory = {'core': lambda x: isinstance(x, list) and len(x) > 0,
						 'core_idx': lambda x: isinstance(x, list) and len(x) > 0,
						 'solver': lambda x: isinstance(x, str)}
		new_optional = {}
		super().__init__()
		self.base_mandatory['method'] = MethodsReconstruction.FASTCORE
		self.add_new_properties(new_mandatory, new_optional)
		self['flux_threshold'] = flux_threshold
		self['core'] = core
		# TODO change this later, this is only for testing
		self['core_idx'] = core
		self['solver'] = solver


class GIMMEProperties(PropertiesReconstruction):
	def __init__(self, exp_vector, objectives, obj_frac=0.9, preprocess=False, flux_threshold=None):
		new_mandatory = {
			'exp_vector': lambda x: isinstance(x, list) and len(x) > 0 or isinstance(x, ndarray),
			'preprocess': lambda x: isinstance(x, bool) or x is None,
			'objectives': lambda x: type(x) in [list, tuple, ndarray]}
		new_optional = {'obj_frac': lambda x: type(x) in [ndarray, list, tuple, float]}
		super().__init__()

		self.add_new_properties(new_mandatory, new_optional)

		self['objectives'] = objectives
		self['exp_vector'] = exp_vector
		self['obj_frac'] = obj_frac if isinstance(obj_frac, ndarray) else array([obj_frac] * len(objectives))
		self['preprocess'] = True if preprocess else False
		self['flux_threshold'] = 1e-4 if flux_threshold is None else flux_threshold


class IMATProperties(PropertiesReconstruction):
	def __init__(self, exp_vector, exp_thresholds, core=None, tolerance=1e-8, epsilon=1):
		new_mandatory = {
			'exp_vector': lambda x: isinstance(x, list) and len(x) > 0 or isinstance(x, ndarray),
			'objectives': lambda x: type(x) in [list, ndarray],
			'exp_thresholds': lambda x: type(x) in (tuple, list, ndarray) and type(x[0]) in [float, int] and type(
				x[1]) in [float, int]
		}
		new_optional = {
			'core': lambda x: type(x) in [ndarray, list, tuple],
			'tolerance': float,
			'epsilon': lambda x: type(x) in [int, float]
		}
		super().__init__()

		self.add_new_properties(new_mandatory, new_optional)

		self['exp_vector'] = exp_vector
		self['exp_thresholds'] = exp_thresholds
		if core:
			self['core'] = core
		if tolerance:
			self['tolerance'] = tolerance
		if epsilon:
			self['epsilon'] = epsilon


class CORDAProperties(PropertiesReconstruction):
	CONSTRAINBY_VAL = 'val'
	CONSTRAINBY_PERC = 'perc'

	def __init__(self, high_conf_rx, medium_conf_rx, neg_conf_rx, pr_to_np=None, constraint=None, constrainby=None,
				 om=None, ntimes=None, nl=None, solver=None, threads=None):
		'''
		:param high_conf_rx: High confidence reactions
		:param medium_conf_rx: Medium confidence reactions
		:param neg_conf_rx: Negative confidence reactions
		:param pr_to_np: Threshold to include NP reactions if PR reactions depend of them
		:param constraint: Constraint value
		:param constrainby: either 'val' (constrain fluxes by value) or 'perc' (constraint by percentage)
		:param om: cost assigned to reactions when calculating dependencies
		:param ntimes: Number of CORSO FBA simulations performed per dependency assessment
		:param nl: Noise added to reaction costs
		'''
		new_mandatory = {k: is_list for k in ['high_conf_rx', 'medium_conf_rx', 'neg_conf_rx']}

		new_optional = {
			# 'met_tests': lambda x: is_list(x) or x is None,
			'pr_to_np': lambda x: isinstance(x, Number),
			'constraint': lambda x: isinstance(x, Number),
			'constrainby': [self.CONSTRAINBY_VAL, self.CONSTRAINBY_PERC],
			'om': lambda x: isinstance(x, Number),
			'ntimes': lambda x: isinstance(x, int) and x > 0,
			'nl': lambda x: isinstance(x, Number) and x >= 0,
			'threads': int
		}

		super().__init__()
		self.add_new_properties(new_mandatory, new_optional)

		vars = [high_conf_rx, medium_conf_rx, neg_conf_rx, pr_to_np, constraint, constrainby, om, ntimes, nl, solver,
				threads]
		defaults = [None, None, None, 2, 1, CORDAProperties.CONSTRAINBY_VAL, 1e4, 5, 1e-2, 'CPLEX', cpu_count() - 1]
		names = ['high_conf_rx', 'medium_conf_rx', 'neg_conf_rx', 'pr_to_np', 'constraint', 'constrainby', 'om',
				 'ntimes', 'nl', 'solver', 'threads']

		for v, k, d in zip(vars, names, defaults):
			self[k] = v if v is not None else d


class tINITProperties(PropertiesReconstruction):
	def __init__(self, reactions_scores, present_metabolites, essential_reactions, production_weight=0.5,
				 allow_excretion=False,
				 no_reverse_loops=False, solver=None):
		# TODO check later if the params input might be necessary to include for the troppo
		new_mandatory = {
			'solver': lambda x: isinstance(x, str)
		}
		new_optional = {
			'reactions_scores': lambda x: is_list_else_empty(x),
			'present_metabolites': lambda x: is_list_else_empty(x),
			'essential_reactions': lambda x: is_list_else_empty(x),
			'production_weight': lambda x: isinstance(x, float),
			'allow_excretion': lambda x: isinstance(x, bool),
			'no_reverse_loops': lambda x: isinstance(x, bool),
		}

		super().__init__()

		self.add_new_properties(new_mandatory, new_optional)

		properties_dict_list = ["reactions_scores", "present_metabolites", "essential_reactions", "production_weight",
								"allow_excretion", "no_reverse_loops", "solver"]
		properties_list = [reactions_scores, if_none_return_list(present_metabolites),
						   if_none_return_list(essential_reactions), if_none_return_list(production_weight),
						   allow_excretion, no_reverse_loops, solver]
		prop_dict = dict(zip(properties_dict_list, properties_list))

		[self.add_if_not_none(*k) for k in prop_dict.items()]


class MBAProperties(PropertiesReconstruction):
	def __init__(self, medium_set, high_set, tolerance, solver):
		new_mandatory = {
			'solver': lambda x: isinstance(x, str),
			'medium_set': lambda x: is_list(x),
			'high_set': lambda x: is_list(x)
		}
		new_optional = {
			'tolerance': lambda x: isinstance(x, float)
		}

		super().__init__()

		self.add_new_properties(new_mandatory, new_optional)

		self['solver'] = solver
		self['medium_set'] = medium_set
		self['high_set'] = high_set
		self['tolerance'] = 1e-8 if tolerance is None else tolerance


class FastCCProperties(PropertiesReconstruction):
	def __init__(self, flux_threshold, method, solver):
		new_mandatory = {
			'flux_threshold': lambda x: isinstance(x, float),
			'method': lambda x: isinstance(x, str),
			'solver': lambda x: isinstance(x, str)
		}
		new_optional = {

		}
		super().__init__()

		self.add_new_properties(new_mandatory, new_optional)

		self['flux_threshold'] = flux_threshold
		if method in ['original', 'nonconvex']:
			self['method'] = method
		else:
			raise Exception('Methods have to be \'original\' or \'nonconvex\'')


if __name__ == '__main__':
	properties = PropertiesReconstruction()
	print(properties.get_mandatory_properties())
	# pro = FastcoreProperties(['a', 'b', 'c'])
	# print(pro.get_mandatory_properties())
	# pro.has_required_properties()
	# # pro['core']

	pro = MBAProperties(['a', 'b'], ['c', 'd'], 1.0, 'cplex')
	print(pro.get_mandatory_properties())
	pro.has_required_properties()
