import abc

from cobamp.utilities.property_management import PropertyDictionary
from cobamp.wrappers.external_wrappers import model_readers


def decode_rx_list(x, rx, mt): return [rx[k] for k in x]


def decode_mt_list(x, rx, mt): return [mt[k] for k in x]


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


class GapfillAlgorithm(object):
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
	"""
	This class is used to define the properties of the reconstruction algorithm.

	Parameters
	----------
	base_mandatory : dict
		Dictionary containing the mandatory properties of the reconstruction algorithm.
		Includes: solver, template_model, omics_type.
	base_optional : dict
		Dictionary containing the optional properties of the reconstruction algorithm.
		Includes: env_conditions, flux_threshold.
	"""

	def __init__(self):
		self.base_mandatory = {'solver': str, 'template_model': model_readers, 'omics_type': 'omics'}
		self.base_optional = {'env_conditions': dict, 'flux_threshold': float}
		super().__init__(self.base_mandatory, self.base_optional)

	@staticmethod
	@abc.abstractmethod
	def from_integrated_scores(scores, **kwargs):
		pass


class GapfillProperties(PropertyDictionary):
	"""
	This class is used to define the properties of the gap-filling algorithms.

	Parameters
	----------
	base_mandatory : dict
		Dictionary containing the mandatory properties of the gap-filling algorithm.
		Includes: solver.
	base_optional : dict
		Dictionary containing the optional properties of the gap-filling algorithm.
		Includes: lsystem_args, avbl_fluxes.
	"""
	__metaclass__ = abc.ABCMeta

	decoder_function = {'lsystem_args': lambda x, r, m: {k: decode_mt_list(v, r, m) for k, v in x.items()},
						'avbl_fluxes': decode_rx_list}

	@abc.abstractmethod
	def __init__(self):
		self.base_mandatory = {'solver': str}
		self.base_optional = {'lsystem_args': dict, 'avbl_fluxes': lambda x: isinstance(x, (list, tuple, set))}
		super().__init__(self.base_mandatory, self.base_optional)
