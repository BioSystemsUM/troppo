from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
import numpy as np

if __name__ == '__main__':
    s_matrix = np.array([[1, -1, 0, 0, -1, 0, -1, 0, 0],
                         [0, 1, -1, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 1, -1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, -1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, -1, 0],
                         [0, 0, 0, 0, 1, 0, 0, 1, -1]])
    M, N = s_matrix.shape

    lower_bounds = np.array([0] * N).astype(float)
    lower_bounds[3] = -10
    upper_bounds = np.array([10] * N).astype(float)

    names = ['R' + str(i + 1) for i in range(N)]
    exp_vector = np.random.random(N)

    properties = GIMMEProperties(exp_vector=exp_vector,
                                 obj_frac=0.8,
                                 objectives=[{8: 1}],
                                 preprocess=True,
                                 flux_threshold=0.8)

    algorithm = GIMME(s_matrix, lower_bounds, upper_bounds, properties)

    solution = algorithm.run()
    print(exp_vector.round(1))
    print(exp_vector > 0.5)
    print(solution)
    print(['R' + str(i + 1) for i, v in enumerate(solution) if v > 0])
