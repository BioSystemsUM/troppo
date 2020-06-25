import abc

from cobamp.utilities.property_management import PropertyDictionary
from cobamp.wrappers.external_wrappers import model_readers

def decode_rx_list(l, rx, mt): return [rx[k] for k in l]
def decode_mt_list(l, rx, mt): return [mt[k] for k in l]


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
	def __init__(self):
		self.base_mandatory = {'solver': str, 'template_model': model_readers,
							   'omics_type': 'omics'}
		self.base_optional = {'env_conditions': dict, 'flux_threshold': float}
		super().__init__(self.base_mandatory, self.base_optional)

	@staticmethod
	@abc.abstractmethod
	def from_integrated_scores(scores, **kwargs):
		pass


class GapfillProperties(PropertyDictionary):
	__metaclass__ = abc.ABCMeta

	decoder_functions = {
		'lsystem_args': lambda x,r,m: {k:decode_mt_list(v, r, m) for k,v in x.items()},
		'avbl_fluxes': decode_rx_list
	}

	@abc.abstractmethod
	def __init__(self):
		self.base_mandatory = {'solver': str}
		self.base_optional = {'lsystem_args': dict, 'avbl_fluxes': lambda x: isinstance(x, (list,tuple,set))}
		super().__init__(self.base_mandatory, self.base_optional)
