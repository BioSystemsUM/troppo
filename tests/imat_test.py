import numpy as np
from troppo.methods.reconstruction.imat import IMAT, IMATProperties

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

    properties = IMATProperties(exp_vector=exp_vector,
                                exp_thresholds=(0.2, 0.5))

    method = IMAT(s_matrix, lower_bounds, upper_bounds, properties)
    print(['R' + str(i + 1) for i in method.run()])
