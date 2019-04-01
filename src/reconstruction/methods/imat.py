import numpy as np
from itertools import chain

from cobamp.core.linear_systems import GenericLinearSystem, VAR_CONTINUOUS, VAR_BINARY
from cobamp.core.optimization import LinearSystemOptimizer


class IMAT():
	def empty_matrix(self, r, c):
		return np.zeros((r,c))

	def __init__(self, S, lb, ub, properties):
		self.S = np.array(S)
		self.lb, self.ub = np.array(lb), np.array(ub)
		self.properties = properties

	def run_imat(self):
		exp_vector = self.properties['exp_vector']
		exp_lb, exp_ub = self.properties['exp_thresholds']
		core = self.properties['core']
		epsilon = self.properties['epsilon']

		high_idx = (np.where(exp_vector >= exp_ub)[0]).astype(int)
		low_idx = (np.where((exp_vector >= 0) & (exp_vector < exp_lb))[0]).astype(int)

		if core:
			high_idx = np.union1d(high_idx, np.array(core))

		lso, lsystem = self.generate_imat_problem(self.S, self.lb, self.ub, high_idx, low_idx, epsilon)

		solution = lso.optimize()
		return solution

	def run(self):
		tol = self.properties['tolerance']
		solution = self.run_imat()
		to_keep = np.where(abs(solution.x())[:self.S.shape[1]] >= tol)[0]

		if solution.status() != 'optimal':
			print('Solution was not optimal')

		return to_keep


	def generate_imat_problem(self, S, lb, ub, high_idx, low_idx, epsilon):
		m,n = S.shape
		nh, nl = len(high_idx), len(low_idx)

		h_ident, l_ident = self.empty_matrix(nh, n), self.empty_matrix(nl, n)
		h_lb, h_ub, l_lb, l_ub = self.empty_matrix(nh,nh), self.empty_matrix(nh,nh), self.empty_matrix(nl,nl), self.empty_matrix(nl,nl)
		h_diag, l_diag = np.diag_indices_from(h_lb), np.diag_indices_from(l_lb)

		if nh > 0:
			h_ident[(np.array(range(nh)), high_idx)] = 1
			h_lb[h_diag] = lb[high_idx] - epsilon
			h_ub[h_diag] = ub[high_idx] + epsilon

		if nl > 0:
			l_ident[(np.array(range(nl)), low_idx)] = 1
			l_lb[l_diag] = lb[low_idx]
			l_ub[l_diag] = ub[low_idx]

		rows = [
			[S, self.empty_matrix(m,nh+nh+nl)],
			[h_ident, h_lb, self.empty_matrix(nh,nh+nl)],
			[h_ident, self.empty_matrix(nh,nh), h_ub, self.empty_matrix(nh,nl)],
			[l_ident, self.empty_matrix(nl,nh*2),l_lb],
			[l_ident, self.empty_matrix(nl,nh*2),l_ub]
		]

		A = np.vstack(list(map(np.hstack,rows)))
		b_lb = [0]*m + list(lb[high_idx]) + [None]*nh + list(lb[low_idx]) + [None]*nl
		b_ub = [0]*m + [None]*nh + list(ub[high_idx]) + [None]*nl + list(ub[low_idx])

		A_lb, A_ub = np.concatenate([lb, np.array([0] * (2 * nh + nl))]), np.concatenate([ub, np.array([1] * (2 * nh + nl))])
		A_vt = [VAR_CONTINUOUS]*n + [VAR_BINARY]*(2*nh+nl)

		## TODO: Move this to optimization on cobamp
		prefix_maker = lambda cd: list([cd[0]+str(i) for i in range(cd[1])])
		A_names = list(chain(*list(map(prefix_maker,[('V',n),('Hpos',nh),('Hneg',nh),('L',nl)]))))

		lsystem = GenericLinearSystem(S=A, var_types=A_vt, lb=A_lb, ub=A_ub, b_lb=b_lb, b_ub=b_ub, var_names=A_names)
		lso = LinearSystemOptimizer(lsystem)

		A_f = np.zeros((A.shape[1]))
		A_f[n:] = 1
		lsystem.set_objective(A_f, False)
		# DEBUG LINES
		#print(high_idx, low_idx)
		#lsystem.write_to_lp('imat.lp')
		return lso, lsystem

