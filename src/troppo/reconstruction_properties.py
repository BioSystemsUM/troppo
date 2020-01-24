from troppo.methods.base import PropertiesReconstruction

from numpy import ndarray, array


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
