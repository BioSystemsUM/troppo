import numpy as np

from cobamp.core.linear_systems import GenericLinearSystem, VAR_CONTINUOUS, VAR_BINARY
from cobamp.core.optimization import LinearSystemOptimizer


class IMAT():
	def empty_matrix(self, r, c):
		return np.zeros((r,c))

	def __init__(self, S, lb, ub, properties):
		self.S = np.array(S)
		self.lb, self.ub = np.array(lb), np.array(ub)
		self.properties = properties

	def run(self):
		exp_vector = self.properties['exp_vector']
		exp_lb, exp_ub = self.properties['exp_thresholds']
		tol = self.properties['tolerance']
		core = self.properties['core']
		epsilon = self.properties['epsilon']

		high_idx = np.where(exp_vector >= exp_ub)[0]
		low_idx = np.where(0 <= exp_vector < exp_lb)[0]

		if core:
			high_idx = np.union1d(high_idx, np.array(core))

		lso, lsystem = generate_imat_problem(self, S, lb, ub, high_idx, low_idx, epsilon)

		solution = lso.optimize()

		to_remove = abs(solution.x()) < tol

		if solution.status() == 'optimal':
			print('Solution was not optimal')

		return to_remove



	def generate_imat_problem(self, S, lb, ub, high_idx, low_idx, epsilon):
		m,n = S.shape
		nh, nl = len(high_idx), len(low_idx)
		h_ident, l_ident = self.empty_matrix(nh, nl), self.empty_matrix(nh, nl)

		h_ident[list(zip(range(nh),high_idx))] = 1
		l_ident[list(zip(range(nl),low_idx))] = 1


		h_lb = np.fromfunction(lambda x,y: lb[high_idx[x]] - epsilon if x==y else 0, (nh, nh))
		h_ub = np.fromfunction(lambda x,y: ub[high_idx[x]] + epsilon if x==y else 0, (nh, nh))

		l_lb = np.fromfunction(lambda x,y: ub[low_idx[x]] if x==y else 0, (nl, nl))
		l_ub = np.fromfunction((lambda x,y: ub[low_idx[x]] if x==y else 0, (nl, nl)))

		rows = [
			[S, self.empty_matrix(m,nh+nh+nl)],
			[h_ident, h_lb, self.empty_matrix(nh,nh+nl)],
			[h_ident, self.empty_matrix(nh,nh), h_ub, self.empty_matrix(nh,nl)],
			[l_ident, self.empty_matrix(nl,nh*2),l_lb],
			[l_ident, self.empty_matrix(nl,nh*2),l_ub]
		]

		filt = lambda lst,idx: list(filter(lambda x: x in idx, lst))

		A = np.vstack(list(map(np.hstack,rows)))
		b_lb = [0]*n + filt(lb, high_idx) + [None]*nh + filt(lb, low_idx) + [None]*nl
		b_ub = [0]*n + [None]*nh + filt(ub, high_idx) + [None]*nl + filt(ub, low_idx)

		A_lb, A_ub = lb + [0]*(2*nh+nl), ub + [1]*(2*nh+nl)
		A_vt = [VAR_CONTINUOUS]*n + [VAR_BINARY]*(2*nh+nl)
		prefix_maker = lambda char, dim: [char+str(i) for i in range(dim)]
		A_names = sum(list(map(prefix_maker,[('V',n),('Hpos',nh),('Hneg',nh),('L',nl)])))

		lsystem = GenericLinearSystem(S=A, var_types=A_vt, lb=A_lb, ub=A_ub, b_lb=b_lb, b_ub=b_ub, var_names=A_names)
		lso = LinearSystemOptimizer(lsystem)

		A_f = np.zeros(len(A.shape[1]))
		A_f[n:] = 1
		lsystem.set_objective(A_f, True)

		return lso, lsystem

