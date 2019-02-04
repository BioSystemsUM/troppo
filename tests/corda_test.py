from reconstruction.methods.corda import CORDA
from reconstruction.reconstruction_properties import CORDAProperties

if __name__ == '__main__':
	import numpy as np

	S = np.array([[1, -1, 0, 0, -1, 0, -1, 0, 0],
				  [0, 1, -1, 0, 0, 0, 0, 0, 0],
				  [0, 1, 0, 1, -1, 0, 0, 0, 0],
				  [0, 0, 0, 0, 0, 1, -1, 0, 0],
				  [0, 0, 0, 0, 0, 0, 1, -1, 0],
				  [0, 0, 0, 0, 1, 0, 0, 1, -1]])
	M, N = S.shape
	lb = np.array([0] * N).astype(float)
	lb[3] = -10
	ub = np.array([10] * N).astype(float)

	names = ['R' + str(i + 1) for i in range(N)]
	plusone = lambda y: [x-1 for x in y]
	properties = CORDAProperties(
		high_conf_rx=plusone([1,2]),
		medium_conf_rx=plusone([3,6]),
		neg_conf_rx=plusone([2,7])
	)

	algorithm = CORDA(S, lb, ub, properties)

	solution = algorithm.run()
	print(solution)
	print(['R'+str(i+1) for i,v in enumerate(solution) if v == 1])
