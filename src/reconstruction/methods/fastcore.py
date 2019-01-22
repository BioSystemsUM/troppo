import numpy as np
from scipy.sparse import csc_matrix, hstack


class FASTcore():
	def __init__(self, S, lb, ub, properties):
		self.S = S
		self.lb, self.ub = lb, ub
		self.properties = properties
		self.n_metabolites, self.n_reactions = self.S.shape
		self.metabolites, self.reactions = np.array([i for i in range(self.n_metabolites)]), np.array(
			[i for i in range(self.n_reactions)])

	def reverse_irreversible_reactions_in_reverse_direction(self):
		'''
		Identifies the irreversible reactions in the reverse direction and returns the S matrix with the signals for the
		metabolites reversed, a vector with the upper bounds as reversed form the lower bounds and a vector with the lower
		bounds as the reverse of the upper ones.
		Returns: S, ub, lb after modifications

		'''
		irrev_reverse_idx = np.where(self.lb <= 0)  # indexes where are the irreversible reverse reactions
		if irrev_reverse_idx.size > 0:  # if they exist
			self.S[:irrev_reverse_idx] = -self.S[:irrev_reverse_idx]
			temp = self.ub[irrev_reverse_idx]
			self.ub[irrev_reverse_idx] = -self.lb[irrev_reverse_idx]
			self.lb[irrev_reverse_idx] = -temp
		return

	def LP7(self, J, basis):
		nJ = J.size # number of irreversible reactions
		# m,n and m2,n2 are the same since the matrix S does not change throughout the algorithm, so we'll be using
		# self.n_metabolites and self.n_reactions

		# objective function
		f = -np.concatenate([np.zeros((self.n_reactions)), np.ones((nJ))])

		# equalities
		# TODO find a way to implement sparse matrix on numpy
		# Aeq = np.concatenate([self.S, np.zeros([self.n_metabolites, nJ])], axis=1)
		Aeq = hstack([self.S, csc_matrix((self.n_reactions, nJ))])
		beq = np.zeros(self.n_metabolites)

		# inequalities
		Ij = csc_matrix((nJ, self.n_reactions))
		# nJ x n_reactions
		Ij[tuple([list(range(nJ)),J])] = -1 # Ij(sub2ind(size(Ij),(1:nj)',J(:))) = -1; from the original code
		Aineq =




		return


	def findSparseMode(self, J, P, singleton, basis=None):
		# epsilon == self.properties['flux_threshold']
		# model, LPProblem == self.S
		if J.size == 0:
			return

		supp = []

		if basis is None:
			basis = []

		if singleton:
			V, basis = self.LP7(J[0], basis)
		else:
			V, basis = self.LP7(J, basis)

	def preprocessing(self):
		self.reverse_irreversible_reactions_in_reverse_direction()
		irreversible_reactions_idx = np.where(self.lb >= 0)

		A = []
		flipped = False
		singleton = False

		J = np.intersect1d(self.properties['core_idx'], irreversible_reactions_idx)  # irreversible reactions in core

		print('J size' + str(len(J)))

		P = np.delete(self.reactions, self.properties['core_idx'], axis=0)  # non-core reactions

		supp, basis = self.findSparseMode(J, P, singleton)

	def fastcore(self):
		pass


if __name__ == '__main__':
	pass
