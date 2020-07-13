from troppo.methods.base import PropertiesReconstruction
from troppo.utilities.list import is_list


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