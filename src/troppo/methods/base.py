import abc

from cobamp.utilities.property_management import PropertyDictionary
from cobamp.wrappers.external_wrappers import model_readers
from troppo.methods_reconstruction import MethodsReconstruction


class ContextSpecificModelReconstructionAlgorithm(object):

	__metaclass__ = abc.ABCMeta
	@abc.abstractmethod
	def __init__(self, S, lb, ub, properties):
		pass

	@abc.abstractmethod
	def run(self):
		pass

	@property
	def properties_class(self):
		return None


class PropertiesReconstruction(PropertyDictionary):
	def __init__(self):
		self.base_mandatory = {'solver': str, 'template_model': model_readers, 'method': MethodsReconstruction,
							   'omics_type': 'omics'}
		self.base_optional = {'env_conditions': dict, 'flux_threshold': float}
		super().__init__(self.base_mandatory, self.base_optional)

	@staticmethod
	@abc.abstractmethod
	def from_integrated_scores(scores, **kwargs):
		pass