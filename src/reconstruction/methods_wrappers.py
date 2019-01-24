import abc
import numpy as np
from cobamp.wrappers.external_wrappers import model_readers
from src.reconstruction.methods.fastcore import FASTcore
from src.reconstruction.reconstruction_properties import FastcoreProperties

map_properties_algorithms = {
	FastcoreProperties : FASTcore
}

class ReconstructionWrapper(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, model):
		self.__model = model
		if model.__module__ in model_readers.keys():
			self.model_reader = model_readers[model.__module__](model)
		else:
			raise TypeError(
				"The `model` instance is not currently supported by cobamp. Currently available readers are: " + str(
					list(model_readers.keys())))
		self.S = self.model_reader.get_stoichiometric_matrix()
		self.lb, self.ub  = [np.array(bounds) for bounds in self.model_reader.get_model_bounds(False, True)]
		# ...
		pass

	def run(self, properties):
		# if properties['core']:
		# 	properties.add_new_properties({'core_idx': lambda x: isinstance(x, list) and len(x) > 0})
		# 	properties['core_idx'] = np.array([self.model_reader.reaction_id_to_index(reaction) for reaction in properties['core']])
		map_properties_algorithms[type(properties)](self.S, self.lb, self.ub, properties)
		# else:
		# 	map_properties_algorithms[type(properties)](self.S, self.lb, self.ub, properties)
		pass

# class FASTcoreWrapper(ReconstructionWrapper):
#
# 	def __init__(self, model, properties):
# 		super.__init__(model)
# 		if isinstance(properties, FastcoreProperties):
# 			self.properties = properties
# 		else:
# 			raise Exception('The properties are not from the FASTcore algorithm')
#
# 	def run(self, proper):
# 		algorithm = FASTcore(self.S, self.lb, self.ub, self.properties)
# 		return [for r in algorithm]



