from troppo.methods.base import PropertiesReconstruction, ContextSpecificModelReconstructionAlgorithm
import numpy as np
import scipy.linalg as linalg
import scipy.sparse as sprs
import scipy.sparse.linalg as sprslinalg

from cobamp.core.linear_systems import fix_backwards_irreversible_reactions, GenericLinearSystem, VAR_CONTINUOUS
from cobamp.core.optimization import LinearSystemOptimizer

EPS = 2^-52

def find_core(S, lb, ub, blocked, weights, solver):
	rev = (lb < 0)
	S = S.copy().astype(float)
	m, n = S.shape
	dense = np.zeros(n,)
	if sum(blocked > 0) > 0:
		dense[blocked > 0] = np.random.normal(size=(sum(blocked > 0),))
	non_null_weights = (abs(weights) > 2^-52)
	k, l = non_null_weights & rev, non_null_weights & np.logical_not(rev)

	objective = np.concatenate([dense, weights[k], weights[l]])
	nk, nl = sum(k), sum(l)
	t1, t2 = np.eye(n), np.eye(nk+nl)

	A = np.vstack([
		np.hstack([S, np.zeros([m, nk+nl])]),
		np.hstack([t1[k,:], t2[rev[non_null_weights],:]]),
		np.hstack([-t1[non_null_weights,:], t2])
	])

	b_lb = np.zeros(m+2*nk+nl,)
	b_ub = b_lb.copy()
	b_ub[m:m+2*nk+nl] = 1e20

	nlb = np.array([-1e20]*n)
	nlb[blocked > 0] = lb[blocked > 0]
	nlb[l] = 0
	if sum(blocked > 0) < 1:
		nlb[l] = 1
	fnlb = np.concatenate([nlb, np.array([-1e20]*(nk+nl))])

	nub = np.array([1e20]*n)
	nub[blocked] = ub[blocked > 0]
	fnub = np.concatenate([nub, np.array([1e20]*(nk+nl))])

	lsys = GenericLinearSystem(A, VAR_CONTINUOUS, fnlb, fnub, b_lb, b_ub,
	                           var_names=['x'+str(i) for i in range(A.shape[1])] ,solver=solver)
	lso = LinearSystemOptimizer(lsys)
	lsys.set_objective(coefficients=objective, minimize=False)
	return lso.optimize()

class SwiftcoreProperties(PropertiesReconstruction):
	def __init__(self, core, weights, flux_threshold=1e-4, solver=None):
		new_mandatory = {'core': lambda x: isinstance(x, list) and len(x) > 0,
		                 'solver': lambda x: isinstance(x, str),
		                 'weights': lambda x: isinstance(x, (list,np.ndarray)) and len(x) > 0}
		new_optional = {'flux_threshold': lambda x: isinstance(x, (int, float))}
		super().__init__()
		# self.base_mandatory['method'] = MethodsReconstruction.FASTCORE
		self.add_new_properties(new_mandatory, new_optional)
		self['flux_threshold'] = flux_threshold
		self['core'] = list(core)
		# TODO change this later, this is only for testing
		self['solver'] = solver
		self['weights'] = weights if weights is not None else np.array([])

	@staticmethod
	def from_integrated_scores(scores, **kwargs):
		return SwiftcoreProperties(core=scores, **dict({k:v for k,v in kwargs.items() if k not in ['core']}))


class SWIFTCORE(ContextSpecificModelReconstructionAlgorithm):
	properties_class = SwiftcoreProperties

	def __init__(self, S, lb, ub, properties):
		super().__init__(S, lb, ub, properties)
		self.S = np.array(S)
		self.model_lb, self.model_ub = np.array(lb), np.array(ub)
		lbnorm, ubnorm  = linalg.norm(self.model_lb, np.inf), linalg.norm(self.model_ub, np.inf)
		self.model_lb, self.model_ub = self.model_lb/lbnorm, self.model_ub/ubnorm
		self.properties = properties
		self.n_metabolites, self.n_reactions = self.S.shape
		self.Sf, self.lbf, self.ubf, fwd_irrev, bwd_irrev = \
			fix_backwards_irreversible_reactions(self.S, self.model_lb, self.model_ub)

	def run(self):
		core = np.array(self.properties['core'])
		tol = self.properties['flux_threshold']
		weights = np.array(self.properties['weights'])
		blocked = np.zeros(len(weights),).astype(bool)
		weights[core] = 0
		flux_sol = find_core(self.Sf, self.lbf, self.ubf, blocked, weights, self.properties['solver'])
		flux = flux_sol.x()[:self.Sf.shape[1]]
		weights[abs(flux) > tol] = 0

		if len(weights) == self.Sf.shape[1]:
			blocked = np.isin(np.arange(self.Sf.shape[1]), core)
			blocked[abs(flux) > tol] = False
		else:
			_, D, V = sprs.linalg.svds(sprs.csc_matrix(self.Sf[:,weights==0]), 10, which='SM')
			wv = linalg.norm(V[:, np.diag(D) < linalg.norm(self.Sf[:,weights == 0], 'fro')*EPS], np.inf ,1)
			blocked[weights == 0] = wv

		while sum(blocked > 0) > 0:
			n_blocked = sum(blocked)
			flux_sol = find_core(self.Sf, self.lbf, self.ubf, blocked, weights, self.properties['solver'])
			flux = flux_sol.x()[:self.Sf.shape[1]]
			weights[abs(flux) > tol] = 0
			blocked[abs(flux) > tol] = 0
			if 2 * sum(blocked) > n_blocked:
				weights = weights / 2

		return np.where(weights == 0)[0]