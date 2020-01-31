from troppo.methods.gimme import GIMME, GIMMEProperties

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
	exp_vector = np.random.random(N)
	exp_vector[[0,2,5,8]] = -1
	properties = GIMMEProperties(
		exp_vector=exp_vector,
		obj_frac=0.8,
		objectives=[np.array([0,0,0,0,0,0,0,0,1])],
		preprocess=True,
		flux_threshold=0.8
	)

	algorithm = GIMME(S, lb, ub, properties)

	solution = algorithm.run()
	print(exp_vector.round(1))
	print(exp_vector > 0.5)
	print(solution)
	print(['R'+str(i+1) for i,v in enumerate(solution) if v > 0])
