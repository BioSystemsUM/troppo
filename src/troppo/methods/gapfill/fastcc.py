# import numpy as np
from numpy import setdiff1d, intersect1d, array, abs, ones, union1d, diag, where, inf, any, vstack, dot
from numpy.linalg import norm
from troppo.utilities.extra_functions_model import ExtraFunctionsModel

from troppo.methods.base import PropertiesReconstruction
from troppo.methods.reconstruction.fastcore import FASTcore


class FastCC(FASTcore):
	# notes: compared to the MATLAB version, modeFlag is always 'on'

	def __init__(self, S, lb, ub, properties):
		super().__init__(S, lb, ub, properties)
		self.epsilon = self.properties['flux_threshold']
		self.method = self.properties['method']
		self.o_S, self.o_lb, self.o_ub = self.S.copy(), self.lb.copy(), self.ub.copy()
		self.f_S, self.f_lb, self.f_ub = None, None, None
		self.I, self.A, self.J = None, None, None
		self.V, self.v, self.Supp = None, None, None
		self.incI, self.irrev_reverse = None, None

	def prepocessing(self):
		extra = ExtraFunctionsModel()
		self.irrev_reverse = where(self.ub <= 0)[0]
		print(self.irrev_reverse)
		self.reverse_irreversible_reactions_in_reverse_direction(self.irrev_reverse)
		self.s_S, self.s_lb, self.s_ub = self.S.copy(), self.lb.copy(), self.ub.copy()
		self.I = where(self.o_lb >= 0)[0]
		self.A = []
		self.J = intersect1d(self.reactions, self.I)
		print(str(self.n_reactions) + ' Total reactions')
		print(str((self.n_reactions - self.I)) + ' Reversible reactions')
		print(str(self.I.size) + ' Irreversible reactions')
		self.V = array([])
		# self.generate_base_LPproblem()
		# self.update_LP7_problem(self.J)
		self.v = self.LP7(self.J, self.epsilon)
		self.Supp = array([i for i, k in self.v.items() if
						   (abs(k) >= 0.99 * self.properties['flux_threshold']) and i <= self.n_reactions - 1])

		self.A = self.Supp
		print(str(self.A.size) + ' Flux consistent reactions, without flipping')

		if self.A.size > 0:
			if self.V.size == 0:
				self.V = array(list(self.v.values()))
			else:
				self.V = vstack([self.V, array(list(self.v.values()))])

		self.incI = setdiff1d(self.J, self.A)  # set of irreversible reaction with inconsistent flux
		if self.incI.size > 0:
			print(str(self.incI.size) + ' Flux inconsistent irreversible reactions, without flipping')

		self.J = setdiff1d(setdiff1d(self.reactions, self.A),
						   self.incI)  # set of reactions with absolute value less than epsilon in the first LP7
		print(str(self.J.size) + ' Flux inconsistent reactions, without flipping')

	def fastcc(self):
		flipped = False
		singleton = False
		JiRev = array([])
		orientation = ones(self.o_S.shape[1]).T
		self.generate_base_LPproblem()
		self.generate_LP3_problem()
		while self.J.size != 0:
			if self.method == 'original':
				if singleton:
					Ji = self.J[0]
					self.v = self.LP3(Ji)
				else:
					Ji = self.J
					self.v = self.LP7(Ji, self.epsilon)
			else:  # TODO implement the fastcc_nonconvex_check_consistency_one_reaction() and fastcc_nonconvex_maximise_card_J
				pass

			# Supp is the set of reactions in v with absolute value greater than epsilon
			self.Supp = array([i for i, k in self.v.items() if
							   (abs(k) >= 0.99 * self.properties['flux_threshold']) and i <= self.n_reactions - 1])

			# A is the set of reactions in self.v with an absolute value greater than epsilon
			self.nA1 = self.A.size
			self.A = union1d(self.A, self.Supp)
			self.nA2 = self.A.size

			# save self.v if new flux consistent reaction found
			if self.nA2 > self.nA1:
				if JiRev.size > 0:
					# l = orientation.size # TODO check if I really need this
					vf = dot(diag(orientation), array(list(self.v.values())))
					self.V = vstack([self.V, vf])

					# sanity check
					print('sanity check')
					if norm(x=dot(self.s_S, vf)) > (self.epsilon / 100):
						print(str(self.epsilon / 100) + ' = epsilon/100')
						print('Should be zero ' + str(norm(x=dot(self.o_S, array(list(self.v.values()))))))
						print('Should be zero ' + str(norm(x=dot(self.s_S, vf))))
						print('May not be zero ' + str(norm(x=dot(self.o_S, array(list(self.v.values()))))))
						print('May not be zero ' + str(norm(x=dot(self.s_S, vf))))
						raise Exception('Flipped flux consistency step failed')
				else:
					self.V = vstack([self.V, array(list(self.v.values()))])
				print(str(self.A.size) + ' Flux consistent reactions')

			# second part -if the set of reactions in V with absolute value less than epsilon has
			# no reactions in common with the set of reactions in V with absolute value
			# greater than epsilon, then flip the sign of the reactions with absolute
			# value less than epsilon because perhaps they are flux consistent in
			# the reverse direction

			if intersect1d(self.J, self.A).size > 0:
				self.J = setdiff1d(self.J, self.A)
				print(str(self.J.size) + ' Flux inconsistent reversible reactions left to flip')
				flipped = False
			else:
				# do not flip the direction of exclusively forward reactions
				JiRev = setdiff1d(Ji, self.I)
				if flipped or JiRev.size == 0:
					flipped = False
					if singleton:
						self.J = setdiff1d(self.J, Ji)
						print(str(self.reactions[Ji]) + ' flux is inconsistent')
					else:
						singleton = True
				else:
					self.reverse_irreversible_reactions_in_reverse_direction_LP3problem(JiRev)
					flipped = True
					orientation[JiRev] = dot(orientation[JiRev], -1)
					print(str(JiRev.size) + ' reversible reactions flipped')

		self.f_S, self.f_lb, self.f_ub = self.S.copy(), self.lb.copy(), self.ub.copy()

		flippedReverseOrientation = ones(self.S.shape[1]).T
		flippedReverseOrientation[self.irrev_reverse] = -1
		self.V = dot(diag(flippedReverseOrientation), self.V.T)

		if norm(dot(self.o_S, self.V), inf) > (self.epsilon / 100):
			print(str((self.epsilon / 100)) + ' = epsilon=100')
			print(str(norm(dot(self.s_S, self.V), inf)) + ' = ||S*V||.')
			print('Flux consistency numerically challenged')
			return self.A, self.f_S, self.f_lb, self.f_ub, self.V
		else:
			print('Flux consistency check finished')
			print((sum(any(
				abs(self.V)) >= 0.99 * self.epsilon)) + ' = Number of flux consistent columns')  # this should fuck up here
			print((norm(dot(self.o_S, self.V)), inf) + ' ||S*V||')
		if self.A.size == self.n_reactions:
			print('Fastcc: the input model is entirely flux consistent')

		return self.A, self.f_S, self.f_lb, self.f_ub, self.V


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