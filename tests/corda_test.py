from troppo.methods.reconstruction.corda import CORDA, CORDAProperties
import numpy as np

if __name__ == '__main__':
    s_matrix = np.array([[1, -1, 0, 0, -1, 0, -1, 0, 0],
                         [0, 1, -1, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 1, -1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, -1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, -1, 0],
                         [0, 0, 0, 0, 1, 0, 0, 1, -1]])
    m, n = s_matrix.shape

    lower_bound = np.array([0] * n).astype(float)
    lower_bound[3] = -10
    upper_bound = np.array([10] * n).astype(float)

    names = ['R' + str(i + 1) for i in range(n)]
    plus_one = lambda y: [x - 1 for x in y]

    properties = CORDAProperties(high_conf_rx=plus_one([1, 2]),
                                 medium_conf_rx=plus_one([3, 6]),
                                 neg_conf_rx=plus_one([2, 7]))

    algorithm = CORDA(s_matrix, lower_bound, upper_bound, properties)

    solution = algorithm.run()
    print(solution)
    print(['R' + str(i + 1) for i, v in enumerate(solution) if v == 1])
