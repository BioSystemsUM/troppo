from cobamp.utilities.property_management import PropertyDictionary
from cobamp.wrappers.external_wrappers import model_readers
from src.reconstruction import MethodsReconstruction


class PropertiesReconstruction(PropertyDictionary):
	def __init__(self):
		self.base_mandatory = {'solver': str, 'template_model': model_readers, 'method': MethodsReconstruction,
							   'omics_type': 'omics'}
		self.base_optional = {'env_conditions': dict, 'flux_threshold': float}
		super().__init__(self.base_mandatory, self.base_optional)


class FastcoreProperties(PropertiesReconstruction):
	def __init__(self, core, flux_threshold=1e-4):
		new_mandatory = {'core': lambda x: isinstance(x, list) and len(x) > 0}
		new_optional = {}
		super().__init__()
		self.base_mandatory['method'] = MethodsReconstruction.FASTCORE
		self.add_new_properties(new_mandatory, new_optional)
		self['flux_threshold'] = flux_threshold
		self['core'] = core


if __name__ == '__main__':
	properties = PropertiesReconstruction()
	print(properties.get_mandatory_properties())
	pro = FastcoreProperties(['a','b','c'])
	print(pro.get_mandatory_properties())
	pro.has_required_properties()
	pro['core']
