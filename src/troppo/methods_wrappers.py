import abc
import numpy as np
from cobamp.wrappers.external_wrappers import model_readers, AbstractObjectReader

from .methods.fastcore import FASTcore
from .methods.gimme import GIMME
from .methods.imat import IMAT
from .methods.corda import CORDA
from .methods.tINIT import tINIT

from .reconstruction_properties import FastcoreProperties, GIMMEProperties, IMATProperties, tINITProperties, \
	CORDAProperties

map_properties_algorithms = {
	FastcoreProperties : FASTcore,
	GIMMEProperties : GIMME,
	IMATProperties : IMAT,
	tINITProperties : tINIT,
	CORDAProperties : CORDA
}

class ReconstructionWrapper(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, model):
		self.__model = model
		if model.__module__ in model_readers.keys():
			self.model_reader = model_readers[model.__module__](model)
		elif isinstance(model, AbstractObjectReader):
			self.model_reader = model
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
		algo = map_properties_algorithms[type(properties)](self.S, self.lb, self.ub, properties)
		return algo.run()
		# else:
		# 	map_properties_algorithms[type(properties)](self.S, self.lb, self.ub, properties)


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