from cobamp.utilities.property_management import PropertyDictionary
from cobamp.wrappers.external_wrappers import model_readers

from numbers import Number
from reconstruction.methods_reconstruction import MethodsReconstruction
from numpy import ndarray, array



class PropertiesReconstruction(PropertyDictionary):
	def __init__(self):
		self.base_mandatory = {'solver': str, 'template_model': model_readers, 'method': MethodsReconstruction,
							   'omics_type': 'omics'}
		self.base_optional = {'env_conditions': dict, 'flux_threshold': float}
		super().__init__(self.base_mandatory, self.base_optional)


class FastcoreProperties(PropertiesReconstruction):
	def __init__(self, core, flux_threshold=1e-4):
		new_mandatory = {'core': lambda x: isinstance(x, list) and len(x) > 0,
						 'core_idx': lambda x: isinstance(x, list) and len(x) > 0}
		new_optional = {}
		super().__init__()
		self.base_mandatory['method'] = MethodsReconstruction.FASTCORE
		self.add_new_properties(new_mandatory, new_optional)
		self['flux_threshold'] = flux_threshold
		self['core'] = core
		# TODO change this later, this is only for testing
		self['core_idx'] = core
	# self['core_idx'] = [model_readers.reaction_id_to_index(reaction) for reaction in core]


class GIMMEProperties(PropertiesReconstruction):
	def __init__(self, exp_vector, objectives, obj_frac=0.9, preprocess=False, flux_threshold=None):
		new_mandatory = {
			'exp_vector': lambda x: isinstance(x, list) and len(x) > 0 or isinstance(x, ndarray),
			'preprocess': lambda x: isinstance(x, bool) or x is None,
			'objectives': lambda x: type(x) in [list, ndarray]}
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
				 om=None, ntimes=None, nl=None):
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
		new_mandatory = {k: is_list for k in ['high_conf_rx','medium_conf_rx','neg_conf_rx']}

		new_optional = {
			#'met_tests': lambda x: is_list(x) or x is None,
			'pr_to_np': lambda x: isinstance(x,Number),
			'constraint': lambda x: isinstance(x,Number),
			'constrainby': [self.CONSTRAINBY_VAL, self.CONSTRAINBY_PERC],
			'om': lambda x: isinstance(x,Number),
			'ntimes': lambda x: isinstance(x, int) and x > 0,
			'nl': lambda x: isinstance(x, Number) and x >= 0
		}

		super().__init__()
		self.add_new_properties(new_mandatory, new_optional)

		vars = [high_conf_rx, medium_conf_rx, neg_conf_rx, pr_to_np, constraint, constrainby, om, ntimes, nl]
		defaults = [None, None, None, 2, 1, CORDAProperties.CONSTRAINBY_VAL, 1e4, 5, 1e-2]
		names = ['high_conf_rx','medium_conf_rx','neg_conf_rx', 'pr_to_np', 'constraint', 'constrainby', 'om', 'ntimes', 'nl']

		for v,k,d in zip(vars,names,defaults):
			self[k] = v if v is not None else d


class tINITProperties(PropertiesReconstruction):
	def __init__(self, reactions_scores, present_metabolites, essential_reactions, production_weight, allow_excretion,
				 no_reverse_loops, params):
		new_mandatory = {
			'reactions_scores': lambda x: is_list(x),
			'present_metabolites': lambda x: is_list(x),
		}
		new_optional = {}

		super.__init__()


		self.add_new_properties(new_mandatory, new_optional)


def is_list(x):
	return type(x) in [list, tuple] and len(x) > 0 or isinstance(x, ndarray) and x.size > 0

def is_number(x):
	return


if __name__ == '__main__':
	properties = PropertiesReconstruction()
	print(properties.get_mandatory_properties())
	pro = FastcoreProperties(['a', 'b', 'c'])
	print(pro.get_mandatory_properties())
	pro.has_required_properties()
	pro['core']
