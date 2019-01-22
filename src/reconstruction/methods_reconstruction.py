from enum import Enum, unique
#from reconstruction.methods import fastcore


@unique
class MethodsReconstruction(Enum):
	MBA = 'mba'
	FASTCORE = 'fastcore'
	INIT = 'init'
	tINIT = 'tinit'
	mCADRE = 'mcadre'
	IMAT = 'imat'
	GIMME = 'gimme'
	EFLUX = 'eflux'

	def describe(self):
		return self.name, self.value

	def get_method(self):
		return self.value


if __name__ == '__main__':
	pass
