from urllib.request import urlretrieve

import numpy as np
from cobra.io import read_sbml_model
from scipy.sparse import csc_matrix, hstack, eye, vstack
from cobamp.core.linear_systems import GenericLinearSystem, VAR_CONTINUOUS
from cobamp.core.optimization import LinearSystemOptimizer
from cobamp.wrappers import COBRAModelObjectReader


class FASTcore():
	def __init__(self, S, lb, ub, properties):
		self.S = S
		self.lb, self.ub = np.array(lb), np.array(ub)
		self.properties = properties
		self.n_metabolites, self.n_reactions = self.S.shape
		self.metabolites, self.reactions = np.array([i for i in range(self.n_metabolites)]), np.array(
			[i for i in range(self.n_reactions)])
		self.counter = 0

	def reverse_irreversible_reactions_in_reverse_direction(self, irrev_reverse_idx):
		'''
		Identifies the irreversible reactions in the reverse direction and returns the S matrix with the signals for the
		metabolites reversed, a vector with the upper bounds as reversed form the lower bounds and a vector with the lower
		bounds as the reverse of the upper ones.
		Returns: S, ub, lb after modifications

		'''
		if irrev_reverse_idx.size > 0:  # if they exist
			self.S[:, irrev_reverse_idx] = -self.S[:, irrev_reverse_idx]
			temp = self.ub[irrev_reverse_idx]
			self.ub[irrev_reverse_idx] = -self.lb[irrev_reverse_idx]
			self.lb[irrev_reverse_idx] = -temp
		return

	def LP7(self, J, basis=None):
		# TODO implement basis (from CPLEX) to improve the speed of the algorithm
		nJ = J.size  # number of irreversible reactions
		# m,n and m2,n2 are the same since the matrix S does not change throughout the algorithm, so we'll be using
		# self.n_metabolites and self.n_reactions

		# objective function
		f = -np.concatenate([np.zeros((self.n_reactions)), np.ones((nJ))])

		# equalities
		Aeq = hstack([self.S, csc_matrix((self.n_metabolites, nJ))])
		beq = np.zeros(self.n_metabolites)

		# inequalities
		Ij = csc_matrix((nJ, self.n_reactions))
		# nJ x n_reactions
		Ij[tuple([list(range(nJ)), J])] = -1  # Ij(sub2ind(size(Ij),(1:nj)',J(:))) = -1; from the original code
		Aineq = hstack([Ij, eye(nJ)])
		bineq = np.zeros((nJ,))

		# bounds
		lb = np.concatenate([self.lb, -np.inf * np.ones((nJ,))])
		ub = np.concatenate([self.ub, np.ones((nJ,)) * self.properties['flux_threshold']])

		A = vstack([Aeq, Aineq])
		b_lb, b_ub = np.concatenate([beq, [None] * nJ]), np.concatenate([beq, bineq])
		# b_lb, b_ub = np.concatenate([beq, -np.inf * np.ones((nJ, ))]), np.concatenate([beq, bineq])

		problem = GenericLinearSystem(A.toarray(), VAR_CONTINUOUS, lb, ub, b_lb, b_ub,
									  ['V' + str(i) for i in range(A.shape[1])])

		lso = LinearSystemOptimizer(problem)
		problem.set_objective(f, True)
		solution = lso.optimize()

		if solution.status() != 'optimal':
			print('Warning, Solution is not optimal')

		if solution:
			return {i:k[1] for i,k in enumerate(solution.var_values().items()) if i<= self.n_reactions-1}
		else:
			return [np.nan] * self.n_reactions

	# return solution  # , problem.model.problem.solution.basis.get_basis()

	def LP9(self, K, P):
		scalingFactor = 1e5

		V = []
		if K.size == 0 or P.size == 0:
			return np.array([])

		nP = P.size
		nK = K.size

		# objective
		f = np.concatenate([np.zeros((self.n_reactions)), np.ones((nP))])

		# equalities
		Aeq = hstack([self.S, csc_matrix((self.n_metabolites, nP))])
		beq = np.zeros(self.n_metabolites)
		# inequalities
		Ip = csc_matrix((nP, self.n_reactions))
		Ip[tuple([list(range(nP)), P])] = 1
		Ik = csc_matrix((nK, self.n_reactions))
		Ik[tuple([list(range(nK)), K])] = 1

		Aineq = vstack([hstack([Ip, -eye(nP)]), hstack([-Ip, -eye(nP)]), hstack([-Ik, csc_matrix((nK, nP))])])
		bineq = np.concatenate(
			[np.zeros((2 * nP,)), -np.ones((nK,)) * self.properties['flux_threshold'] * scalingFactor])
		# bounds
		lb = np.concatenate([self.lb, np.zeros((nP,))]) * scalingFactor
		ub = np.concatenate([self.ub, np.array(list(map(max, zip(np.abs(self.lb[P]), np.abs(self.ub[P])))))]) * scalingFactor

		A = vstack([Aeq, Aineq])
		b_lb, b_ub = np.concatenate([beq, [None] * (2 * nP + nK)]), np.concatenate([beq, bineq])
		# b_lb, b_ub = np.concatenate([beq, -np.inf * np.ones((nJ, ))]), np.concatenate([beq, bineq])

		problem = GenericLinearSystem(A.toarray(), VAR_CONTINUOUS, lb, ub, b_lb, b_ub,
									  ['V' + str(i) for i in range(A.shape[1])])

		lso = LinearSystemOptimizer(problem)
		problem.set_objective(f, True)
		# problem.write_to_lp('Test_LP9_'+str(self.counter))
		# self.counter+=1
		solution = lso.optimize()

		if solution.status() != 'optimal':
			print('Warning, Solution is not optimal')

		if solution:
			return {i:k[1] for i,k in enumerate(solution.var_values().items()) if i<= self.n_reactions-1}
		else:
			return [np.nan] * self.n_reactions

	def findSparseMode(self, J, P, singleton, basis=None):
		# epsilon == self.properties['flux_threshold']
		# model, LPProblem == self.S
		if J.size == 0:
			return np.array([], dtype=int)

		supp = []

		if basis is None:
			basis = []

		if singleton:
			# V, basis = self.LP7(J[0], basis)
			V = self.LP7(J[0], basis)
		else:
			# V, basis = self.LP7(J, basis)
			V = self.LP7(J, basis)

		# K = np.array([rJ for rJ in J if V['V' + str(rJ)] >= 0.99 * self.properties['flux_threshold']])
		K = np.array([rJ for rJ in J if V[rJ] >= 0.99 * self.properties['flux_threshold']])

		print('done LP7')
		if K.size > 0:
			V = self.LP9(K, P)
			Supp = np.array([i for i, k in V.items() if
							 (np.abs(k) >= 0.99 * self.properties['flux_threshold']) and i <= self.n_reactions-1])
			# Supp = np.array([k for k, v in V.items() if v >= 0.99 * self.properties['flux_threshold']])
			print('done LP9')
			return Supp
		else:
			return np.array([])

	def preprocessing(self):
		irreversible_reactions_idx_to_change = np.where(self.ub <= 0)[0]
		self.reverse_irreversible_reactions_in_reverse_direction(irreversible_reactions_idx_to_change)
		irreversible_reactions_idx = np.where(self.lb >= 0)[0]

		singleton = False

		J = np.intersect1d(self.properties['core_idx'], irreversible_reactions_idx)  # irreversible reactions in core

		print('J size' + str(len(J)))
		print(J)

		P = np.setdiff1d(self.reactions, self.properties['core_idx'])  # non-core reactions

		# supp, basis = self.findSparseMode(J, P, singleton)
		supp = self.findSparseMode(J, P, singleton)

		if np.setdiff1d(J, supp).size > 0:
			print('Inconsistent irreversible core reactions')

		return np.setdiff1d(self.properties['core_idx'], supp), supp, P, irreversible_reactions_idx
		# return np.setdiff1d(J, supp), supp, P, irreversible_reactions_idx #TODO this is a test

	def fastcore(self):
		flipped = False
		singleton = False
		J, A, P, irreversible_reactions_idx = self.preprocessing()

		while J.size > 0:
			# print(J)
			P = np.setdiff1d(P, A)

			supp = self.findSparseMode(J, P, singleton)

			# print(A, supp)

			A = np.union1d(A, supp)

			if np.intersect1d(J, A) != np.array([]):
				J = np.setdiff1d(J, A)
				flipped = False
			else:
				# Sara modification
				# if flipped:
				# 	flipped = False
				# 	singleton = True
				#
				if singleton:
					JiRev = np.setdiff1d(J[0], irreversible_reactions_idx)
				else:
					JiRev = np.setdiff1d(J, irreversible_reactions_idx)
				if flipped or JiRev == np.array([]):
					if singleton:
						print('Error: Global network is not consistent')
						return sorted(A)
					else:
						flipped = False
						singleton = True
				else:
					self.reverse_irreversible_reactions_in_reverse_direction(JiRev)
					flipped = True
					print('Flipped')

		# core_reactions_bool = np.array([False] * (self.n_reactions+1))
		# core_reactions_bool[A] = 1
		#
		# tissue_reactions = np.delete(self.reactions, A)

		return sorted(A)


if __name__ == '__main__':
	from reconstruction.reconstruction_properties import FastcoreProperties

	path, content = urlretrieve('http://bigg.ucsd.edu/static/models/iAF1260.xml')

	model = read_sbml_model(path)
	cobamp_model = COBRAModelObjectReader(model)

	cobamp_model.get_irreversibilities(True)

	S = cobamp_model.get_stoichiometric_matrix()
	lb, ub = cobamp_model.get_model_bounds(False, True)
	rx_names = cobamp_model.get_metabolite_and_reactions_ids()[0]

	# S = np.array([[1, -1, 0, 0, -1, 0, -1, 0, 0],
	# 			  [0, 1, -1, 0, 0, 0, 0, 0, 0],
	# 			  [0, 1, 0, 1, -1, 0, 0, 0, 0],
	# 			  [0, 0, 0, 0, 0, 1, -1, 0, 0],
	# 			  [0, 0, 0, 0, 0, 0, 1, -1, 0],
	# 			  [0, 0, 0, 0, 1, 0, 0, 1, -1]])
	# M, N = S.shape
	# lb = np.array([0] * N).astype(float)
	# lb[3] = -10
	# ub = np.array([10] * N).astype(float)
	#
	# names = ['R' + str(i + 1) for i in range(N)]

	f = FASTcore(S, lb, ub, FastcoreProperties(core=[142,133,573, 354, 656, 128,1200]))
	# f = FASTcore(S, lb, ub, FastcoreProperties(core=[8]))
	tissue_reactions = f.fastcore()
	print([i + 1 for i in tissue_reactions])
