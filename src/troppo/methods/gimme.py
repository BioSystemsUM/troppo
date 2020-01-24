import numpy as np
from numpy.core._multiarray_umath import ndarray, array

from cobamp.core.models import GIMMEModel, ConstraintBasedModel
from troppo.methods.base import ContextSpecificModelReconstructionAlgorithm, PropertiesReconstruction

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

	@staticmethod
	def from_integrated_scores(scores, **kwargs):
		return GIMMEProperties(exp_vector=scores, **{k:v for k,v in kwargs.items() if 'exp_vector' not in k})

class GIMME(ContextSpecificModelReconstructionAlgorithm):
	properties_class = GIMMEProperties

	def __init__(self, S, lb, ub, properties):
		super().__init__(S, lb, ub, properties)
		self.S = np.array(S)
		self.lb, self.ub = np.array(lb), np.array(ub)
		self.properties = properties
		self.model = GIMMEModel
		self.sol = None
		metabolite_names = ['M'+str(i) for i in range(S.shape[0])]
		reaction_names = ['R'+str(i) for i in range(S.shape[1])]
		cbm = ConstraintBasedModel(S, list(zip(lb,ub)), reaction_names=reaction_names, metabolite_names=metabolite_names)
		self.gm = GIMMEModel(cbm, self.properties['solver'])


	def run(self):
		sol = self.gm.optimize_gimme(
			exp_vector=self.properties['exp_vector'],
			objectives=self.properties['objectives'],
			obj_frac=self.properties['obj_frac'],
			flux_thres=self.properties['flux_threshold']
		)
		self.sol = sol
		return sol.get_reaction_activity(self.properties['flux_threshold'])

