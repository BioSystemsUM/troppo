import numpy as np
from itertools import chain

from numpy.core._multiarray_umath import ndarray

from cobamp.core.linear_systems import GenericLinearSystem, VAR_CONTINUOUS, VAR_BINARY
from cobamp.core.optimization import LinearSystemOptimizer
from troppo.methods.base import ContextSpecificModelReconstructionAlgorithm, PropertiesReconstruction


class IMATProperties(PropertiesReconstruction):
    """
    Properties for IMAT algorithm

    Parameters
    ----------
    exp_vector : np.ndarray or list
        The vector of expression values
    exp_thresholds : tuple
        The thresholds for the expression values
    core : list, optional
        The core reactions, by default None
    tolerance : float, optional
        The tolerance, by default 1e-8
    epsilon : float, optional
        The epsilon, by default 1
    """
    def __init__(self, exp_vector: np.ndarray or list, exp_thresholds: tuple or list or ndarray,
                 core: ndarray or list or tuple = None, tolerance: float = 1e-8, epsilon: int or float = 1):
        new_mandatory = {
            'exp_vector': lambda x: isinstance(x, list) and len(x) > 0 or isinstance(x, ndarray),
            'exp_thresholds': lambda x: type(x) in (tuple, list, ndarray) and type(x[0]) in [float, int] and type(
                x[1]) in [float, int]
        }
        new_optional = {
            'core': lambda x: type(x) in [ndarray, list, tuple],
            'tolerance': float,
            'epsilon': lambda x: type(x) in [int, float]
        }
        super().__init__()

        self.add_new_properties(new_mandatory, new_optional)

        self['exp_vector'] = exp_vector
        self['exp_thresholds'] = exp_thresholds
        if core:
            self['core'] = core
        if tolerance:
            self['tolerance'] = tolerance
        if epsilon:
            self['epsilon'] = epsilon

    @staticmethod
    def from_integrated_scores(scores: list, **kwargs) -> 'IMATProperties':
        """
        Create IMAT properties from integrated scores

        Parameters
        ----------
        scores: list
            The list of integrated scores
        kwargs: dict
            Additional arguments

        Returns
        -------
        IMATProperties: The IMAT properties

        """
        return IMATProperties(exp_vector=scores, **kwargs)


class IMAT(ContextSpecificModelReconstructionAlgorithm):
    """
    IMAT algorithm

    Parameters
    ----------
    S : ndarray
        The stoichiometric matrix
    lb : ndarray
        The lower bounds
    ub : ndarray
        The upper bounds
    properties : IMATProperties
        The IMAT properties

    Attributes
    ----------
    sol : LinearSystemOptimizer
        The solution
    """

    properties_class = IMATProperties

    @staticmethod
    def empty_matrix(r, c):
        return np.zeros((r, c))

    def __init__(self, S: ndarray, lb: ndarray, ub: ndarray, properties: IMATProperties):
        super().__init__(S, lb, ub, properties)
        self.S = np.array(S)
        self.lb, self.ub = np.array(lb), np.array(ub)
        self.properties = properties
        self.sol = None

    def run_imat(self):
        """
        Run the IMAT algorithm

        Returns
        -------
        solution: Solution
            Instance with the solution to the IMAT problem
        """
        exp_vector = self.properties['exp_vector']
        exp_lb, exp_ub = self.properties['exp_thresholds']
        core = self.properties['core']
        epsilon = self.properties['epsilon']

        high_idx = (np.where(np.array(exp_vector) >= exp_ub)[0]).astype(int)
        low_idx = (np.where((np.array(exp_vector) >= 0) & (np.array(exp_vector) < exp_lb))[0]).astype(int)

        if core:
            high_idx = np.union1d(high_idx, np.array(core))

        lso, lsystem = self.generate_imat_problem(self.S, self.lb, self.ub, high_idx, low_idx, epsilon)

        solution = lso.optimize()
        return solution

    def run(self):
        """
        Run the algorithm and return a list of reactions to keep in the final model.

        Returns
        -------
        list: The list of reactions to keep in the final model

        """
        tol = self.properties['tolerance']
        solution = self.run_imat()
        to_keep = np.where(abs(solution.x())[:self.S.shape[1]] >= tol)[0]

        self.sol = solution

        if solution.status() != 'optimal':
            print('Solution was not optimal')

        return to_keep

    def generate_imat_problem(self, S, lb, ub, high_idx, low_idx, epsilon):
        """
        Generate the IMAT problem

        Parameters
        ----------
        S: ndarray
            The stoichiometric matrix
        lb: ndarray
            The lower bounds
        ub: ndarray
            The upper bounds
        high_idx: ndarray
            The high index
        low_idx: ndarray
            The low index
        epsilon: float or int
            The epsilon value

        Returns
        -------
        lso: LinearSystemOptimizer
            The linear system optimizer
        lsystem: GenericLinearSystem
            The linear system
        """
        m, n = S.shape
        nh, nl = len(high_idx), len(low_idx)

        h_ident, l_ident = self.empty_matrix(nh, n), self.empty_matrix(nl, n)
        h_lb, h_ub, l_lb, l_ub = self.empty_matrix(nh, nh), self.empty_matrix(nh, nh), \
            self.empty_matrix(nl, nl), self.empty_matrix(nl, nl)
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
            [S, self.empty_matrix(m, nh + nh + nl)],
            [h_ident, h_lb, self.empty_matrix(nh, nh + nl)],
            [h_ident, self.empty_matrix(nh, nh), h_ub, self.empty_matrix(nh, nl)],
            [l_ident, self.empty_matrix(nl, nh * 2), l_lb],
            [l_ident, self.empty_matrix(nl, nh * 2), l_ub]
        ]

        A = np.vstack(list(map(np.hstack, rows)))
        b_lb = [0] * m + list(lb[high_idx]) + [None] * nh + list(lb[low_idx]) + [None] * nl
        b_ub = [0] * m + [None] * nh + list(ub[high_idx]) + [None] * nl + list(ub[low_idx])

        A_lb, A_ub = np.concatenate([lb, np.array([0] * (2 * nh + nl))]), np.concatenate(
            [ub, np.array([1] * (2 * nh + nl))])
        A_vt = [VAR_CONTINUOUS] * n + [VAR_BINARY] * (2 * nh + nl)

        ## TODO: Move this to optimization on cobamp
        prefix_maker = lambda cd: list([cd[0] + str(i) for i in range(cd[1])])
        A_names = list(chain(*list(map(prefix_maker, [('V', n), ('Hpos', nh), ('Hneg', nh), ('L', nl)]))))

        lsystem = GenericLinearSystem(S=A, var_types=A_vt, lb=A_lb, ub=A_ub, b_lb=b_lb, b_ub=b_ub, var_names=A_names)
        lso = LinearSystemOptimizer(lsystem)

        A_f = np.zeros((A.shape[1]))
        A_f[n:] = 1
        lsystem.set_objective(A_f, False)
        # DEBUG LINES
        # print(high_idx, low_idx)
        # lsystem.write_to_lp('imat.lp')
        return lso, lsystem
