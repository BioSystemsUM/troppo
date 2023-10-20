import numpy as np
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties

if __name__ == '__main__':
    s_matrix = np.array([[1, -1, 0, 0, -1, 0, -1, 0, 0],
                         [0, 1, -1, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 1, -1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, -1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, -1, 0],
                         [0, 0, 0, 0, 1, 0, 0, 1, -1]])
    M, N = s_matrix.shape

    lower_bounds = np.array([0] * N).astype(float)
    lower_bounds[3] = -1000
    upper_bounds = np.array([1000] * N).astype(float)

    names = ['R' + str(i + 1) for i in range(N)]

    asd = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    asd[0] = 1

    asd = [0] * 9
    asd[8] = 1

    tinit_algorithm = tINIT(s_matrix, lower_bounds, upper_bounds, tINITProperties(reactions_scores=asd,
                                                                                  present_metabolites=[],
                                                                                  essential_reactions=[8],
                                                                                  production_weight=0.0,
                                                                                  allow_excretion=False,
                                                                                  no_reverse_loops=True))

    tinit_algorithm.preprocessing()
    tinit_algorithm.build_problem()

    result = tinit_algorithm.solve_problem()
    result.sort()
    print(result + 1)
