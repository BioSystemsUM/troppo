if __name__ == '__main__':
	import numpy as np
	from numpy import int_
	from reconstruction.methods.tINIT import tINIT
	from reconstruction.reconstruction_properties import tINITProperties

	S = np.array([[1, -1, 0, 0, -1, 0, -1, 0, 0],
				  [0, 1, -1, 0, 0, 0, 0, 0, 0],
				  [0, 1, 0, 1, -1, 0, 0, 0, 0],
				  [0, 0, 0, 0, 0, 1, -1, 0, 0],
				  [0, 0, 0, 0, 0, 0, 1, -1, 0],
				  [0, 0, 0, 0, 1, 0, 0, 1, -1]])
	M, N = S.shape
	lb = np.array([0] * N).astype(float)
	lb[3] = -1000
	# lb[4] = -1000
	# lb[5] = -1000
	ub = np.array([1000] * N).astype(float)

	names = ['R' + str(i + 1) for i in range(N)]

	asd = [0, 0, 0, 0, 0, 0, 0, 0, 0]
	asd[0] = 1

	asd = [0] * 9
	asd[8] = 1
	# asd[]

	t = tINIT(S, lb, ub, tINITProperties(reactions_scores=asd, present_metabolites=[], essential_reactions=[8],
										 production_weight=0.0, allow_excretion=False, no_reverse_loops=True))
	t.preprocessing()
	t.build_problem()
	res = t.solve_problem()
	res.sort()
	print(res+1)
