from multiprocessing import cpu_count
from numbers import Number

from numpy import array, zeros, sqrt, vstack, where, floor, random, unique, \
	apply_along_axis, nan, logical_or

from cobamp.core.models import CORSOModel
from cobamp.core.models import ConstraintBasedModel
from time import time

from troppo.methods.base import ContextSpecificModelReconstructionAlgorithm, PropertiesReconstruction

from pathos.multiprocessing import cpu_count
from pathos.pools import _ProcessPool

from troppo.reconstruction_properties import is_list


def _init_corda_worker(corso_fba, constraint, constrainby, costfx, costbase, ntimes, eps, lb):
	global _corso_fba, _constraint, _constrainby, _costfx, _costbase, _ntimes, _eps, _lb
	_corso_fba = corso_fba
	_constraint = constraint
	_constrainby = constrainby
	_costfx = costfx
	_costbase = costbase
	_ntimes = ntimes
	_eps = eps
	_lb = lb

def _corda_dependent_reactions_iteration(rx):
	global _corso_fba, _constraint, _constrainby, _costfx, _costbase
	global _ntimes, _forward, _eps, _lb
	dependent, to_delete = _corda_dependent_rxs_per_sense(rx, _constraint, True)

	if _lb[rx] < 0:
		bkw_dep, to_del_bkw = _corda_dependent_rxs_per_sense(rx, -_constraint, False)

		dependent = dependent | bkw_dep
		to_delete = to_del_bkw & to_delete

	return rx, (dependent, to_delete)

def _corda_dependent_rxs_per_sense(rx, constraint, forward):
	global _corso_fba, _constrainby, _costfx, _costbase
	of_dict = {rx: 1}
	cost = _costbase + _costfx()
	flux, corso_sol = _corso_fba.optimize_corso(cost, of_dict, not forward, constraint, _constrainby, eps=_eps)

	dependent = abs(corso_sol.x()) > _eps
	to_del = not dependent.any()
	if not to_del:
		for i in range(_ntimes - 1):
			cost = _costbase + _costfx()
			flux, corso_sol = _corso_fba.optimize_corso(cost, of_dict, not forward, constraint, _constrainby, eps=_eps, flux1=flux)
			dependent = (abs(corso_sol.x()) > _eps) | dependent
	else:
		dependent = zeros(dependent.shape).astype(bool)
	return dependent, to_del

class CORDAProperties(PropertiesReconstruction):
	CONSTRAINBY_VAL = 'val'
	CONSTRAINBY_PERC = 'perc'

	def __init__(self, high_conf_rx, medium_conf_rx, neg_conf_rx, pr_to_np=None, constraint=None, constrainby=None,
				 om=None, ntimes=None, nl=None, solver=None, threads=None):
		'''
		:param high_conf_rx: High confidence reactions
		:param medium_conf_rx: Medium confidence reactions
		:param neg_conf_rx: Negative confidence reactions
		:param pr_to_np: Threshold to include NP reactions if PR reactions depend of them
		:param constraint: Constraint value
		:param constrainby: either 'val' (constrain fluxes by value) or 'perc' (constraint by percentage)
		:param om: cost assigned to reactions when calculating dependencies
		:param ntimes: Number of CORSO FBA simulations performed per dependency assessment
		:param nl: Noise added to reaction costs
		'''
		new_mandatory = {k: is_list for k in ['high_conf_rx', 'medium_conf_rx', 'neg_conf_rx']}

		new_optional = {
			# 'met_tests': lambda x: is_list(x) or x is None,
			'pr_to_np': lambda x: isinstance(x, Number),
			'constraint': lambda x: isinstance(x, Number),
			'constrainby': [self.CONSTRAINBY_VAL, self.CONSTRAINBY_PERC],
			'om': lambda x: isinstance(x, Number),
			'ntimes': lambda x: isinstance(x, int) and x > 0,
			'nl': lambda x: isinstance(x, Number) and x >= 0,
			'threads': int
		}

		super().__init__()
		self.add_new_properties(new_mandatory, new_optional)

		vars = [high_conf_rx, medium_conf_rx, neg_conf_rx, pr_to_np, constraint, constrainby, om, ntimes, nl, solver,
				threads]
		defaults = [None, None, None, 2, 1, CORDAProperties.CONSTRAINBY_VAL, 1e4, 5, 1e-2, 'CPLEX', cpu_count() - 1]
		names = ['high_conf_rx', 'medium_conf_rx', 'neg_conf_rx', 'pr_to_np', 'constraint', 'constrainby', 'om',
				 'ntimes', 'nl', 'solver', 'threads']

		for v, k, d in zip(vars, names, defaults):
			self[k] = v if v is not None else d

	@staticmethod
	def from_integrated_scores(scores, **kwargs):
		hi, med, neg = scores
		return CORDAProperties(hi, med, neg, **kwargs)


class CORDA(ContextSpecificModelReconstructionAlgorithm):
	properties_class = CORDAProperties

	def costfx_factory(self, nl, om, costbase):
		return lambda: nl * floor(om * random.rand(len(costbase),)) / om

	def __init__(self, S, lb, ub, properties):
		super().__init__(S, lb, ub, properties)
		self.S = array(S)
		self.lb, self.ub = array(lb), array(ub)
		self.properties = properties

		rx_names, mt_names = ['V' + str(i) for i in range(S.shape[1])], ['M' + str(i) for i in range(S.shape[0])]
		cbmodel = ConstraintBasedModel(S, list(zip(lb, ub)), reaction_names=rx_names,
									   metabolite_names=mt_names, optimizer=True, solver=properties['solver'])
		self._m, self._n = self.S.shape
		self.corso_fba = CORSOModel(cbmodel, solver=properties['solver'])

	def run(self):
		rx_cat = zeros(self._n,)
		rx_cat[self.properties['high_conf_rx']] = 1
		rx_cat[self.properties['medium_conf_rx']] = 2
		rx_cat[self.properties['neg_conf_rx']] = 3

		constraint = self.properties['constraint']
		constrainby = self.properties['constrainby']
		nl = self.properties['nl']
		ntimes = self.properties['ntimes']
		om = self.properties['om']
		pr_to_np = self.properties['pr_to_np']
		threads = self.properties['threads']

		return self.run_corda(rx_cat, constraint, constrainby, nl, ntimes, om, pr_to_np, threads)

	def run_corda(self, rx_cat, constraint, constrainby, nl, ntimes, om=1e4, pr_to_np=2, threads=cpu_count()-1):
		def _corda_find_all_dependencies(reaction_list):
			res_map = {r:i for i,r in enumerate(reaction_list)}
			true_threads = min((len(reaction_list)//2)+1, threads)
			result = [None] * len(reaction_list)
			rx_per_job = len(reaction_list) // threads
			pool = _ProcessPool(
				processes=true_threads,
				initializer=_init_corda_worker,
				initargs=(self.corso_fba, constraint, constrainby, costfx, costbase, ntimes, 1e-6, self.lb)
			)
			for i, value in pool.imap_unordered(_corda_dependent_reactions_iteration, reaction_list,
												chunksize=rx_per_job):
				result[res_map[i]] = value

			pool.close()
			pool.join()

			return result

		import pandas as pd

		rx_cat = array(rx_cat)
		costbase = zeros(self._n, )

		print(pd.Series(rx_cat).value_counts())

		costbase[rx_cat == 2] = sqrt(om)
		costbase[rx_cat == 3] = om

		costfx = self.costfx_factory(nl, om, costbase)

		def nested_dependent_rxs(rx):
			return self.find_dependent_reactions(rx, constraint, constrainby, costfx, costbase, ntimes, eps=1e-6)

		HC_reactions = where(rx_cat == 1)[0]

		print('Step 1 started')
		s1t = time()


			# print(sum(dep),pd.Series(rx_cat).value_counts())
		if (threads and threads > 1) and (len(HC_reactions)/2 > threads):
			res1 = _corda_find_all_dependencies(HC_reactions)
		else:
			res1 = list(map(nested_dependent_rxs, HC_reactions))

		s1_deps, s1_block = list(zip(*res1))
		rx_cat[logical_or.reduce(s1_deps)] = 1
		if sum(s1_block) > 0:
			rx_cat[HC_reactions[s1_block]] = -1

		print('\t- Finished in ',str(time()-s1t),'seconds')
		print(pd.Series(rx_cat).value_counts())
		# for rx in HC_reactions:
		# 	dep, to_del = self.find_dependent_reactions(rx, constraint, constrainby, costfx, costbase, ntimes, eps=1e-6)
		# 	rx_cat[dep] = 1
		# 	if to_del:
		# 		rx_cat[rx] = -1

		self.block_reactions_from_idxs(rx_cat)


		costbase = zeros(self._n, )
		costbase[rx_cat == 3] = om

		PR_reactions = where(rx_cat == 2)[0]
		NP_reactions = where(rx_cat == 3)[0]

		costfx = self.costfx_factory(nl, om, costbase)

		print('Step 2 started')
		s1t = time()
		PR_NP = {}

		# def _step_two(rx):
		# 	dep, to_del = self.find_dependent_reactions(rx, constraint, constrainby, costfx, costbase, ntimes, eps=1e-6)
		# 	PR_NP[rx] = dep[NP_reactions]
		# 	if to_del:
		# 		rx_cat[rx] = -1

		if (threads and threads > 1) and (len(PR_reactions)/2 > threads):
			res2 = _corda_find_all_dependencies(PR_reactions)
		else:
			res2 = list(map(nested_dependent_rxs, PR_reactions))

		s2_deps, s2_block = list(zip(*res2))

		for rx, dep in zip(PR_reactions, s2_deps):
			PR_NP[rx] = dep[NP_reactions]

		if sum(s2_block) > 0:
			rx_cat[PR_reactions[s2_block]] = -1

		print('\t- Finished in ',str(time()-s1t),'seconds')
		print(pd.Series(rx_cat).value_counts())

		self.block_reactions_from_idxs(rx_cat)

		PR_NP = [PR_NP[k] for k in sorted(PR_NP.keys()) if rx_cat[k] != -1]

		if len(PR_NP) > 0:
			PR_NP = vstack(PR_NP)
			NP_occurrence = apply_along_axis(sum, 0, PR_NP)
			np_to_pr_idx = NP_reactions[NP_occurrence > pr_to_np]

			if len(np_to_pr_idx) > 0:
				rx_cat[np_to_pr_idx] = 2
				PR_NP = PR_NP[:, NP_occurrence > pr_to_np]
				if PR_NP.shape[1] > 0:
					PR_NP = vstack([zeros((len(np_to_pr_idx), PR_NP.shape[1]))])
				else:
					PR_NP = array([])

		# 2.2
		PR_reactions = where(rx_cat == 2)[0]
		NP_reactions = where(rx_cat == 3)[0]

		PR_to_check_l8r = []
		rx_cat[NP_reactions] = -1
		self.block_reactions_from_idxs(rx_cat)

		res2 = []
		for i, rx in enumerate(PR_reactions):
			to_del = self.check_if_blocked(rx)
			if to_del:
				# PR_to_check_l8r.append(rx)
				# if len(PR_NP) > 0 and len(np_to_pr_idx) > 0:
				# 	np_from_rx = where(PR_NP[i, :] > 0)[0]
				# 	if len(np_from_rx) == 0:
				# 		print('Undefined')
				# 	else:
				# 		for kn in np_from_rx:
				# 			res2.append(kn)
				rx_cat[rx] = -1
		res2 = unique(sorted(array(res2)))

		PR_reactions = where(rx_cat == 2)[0]
		rx_cat[PR_reactions] = 1
		# rescued = set(sorted(PR_to_check_l8r + res2.tolist()))

		ES_reactions = rx_cat == 1
		OT_reactions = rx_cat == 0

		#to_block = where(ES_reactions | OT_reactions)[0]

		#rx_cat[to_block] = -1 # TODO: check this
		self.block_reactions_from_idxs(rx_cat)

		costbase = zeros(self._n, )
		costbase[OT_reactions] = om

		print('Step 3 started')
		s1t = time()

		ES_OT = {}


		# def _step_three(rx):
		# 	dep, to_del = self.find_dependent_reactions(rx, constraint, constrainby, costfx, costbase, ntimes, 1e-6)
		# 	ES_OT[rx] = dep[OT_reactions]
		#
		ES_temp = where(ES_reactions)[0]
		if (threads and threads > 1) and (len(ES_temp)/2 > threads):
			res3 = _corda_find_all_dependencies(ES_temp)
		else:
			res3 = list(map(nested_dependent_rxs, ES_temp))



		s3_deps, _ = list(zip(*res3))
		for rx, dep in zip(where(ES_reactions)[0], s3_deps):
			ES_OT[rx] = dep[OT_reactions]

		ES_OT = [ES_OT[k] for k in sorted(ES_OT.keys()) if rx_cat[k] != -1]

		print('\t- Finished in ',str(time()-s1t),'seconds')
		print(pd.Series(rx_cat).value_counts())


		OT_reaction_ids = where(OT_reactions)[0]
		if ES_OT:
			ES_OT = vstack(ES_OT)
			if ES_OT.shape[1] > 0:
				OT_occurrence = apply_along_axis(sum, 0, ES_OT)
				ot_to_es_idx = OT_reaction_ids[OT_occurrence != 0]
				rx_cat[ot_to_es_idx] = 1

		return rx_cat

	def block_reactions_from_idxs(self, rxcat):
		block_corso = lambda rx: self.corso_fba.set_reaction_bounds(rx, lb=0, ub=0)
		block_cbmodel = lambda rx: self.corso_fba.cbmodel.set_reaction_bounds(rx, lb=0, ub=0)
		to_block = where(rxcat == -1)[0]
		for rx in to_block:
			self.do_function_for_reactions_on_both_models(rx, block_cbmodel, block_corso)

	def do_function_for_reactions_on_both_models(self, reaction, mfunction, corsofunction):
		if isinstance(self.corso_fba.mapping[reaction], int):
			corsofunction(reaction)
		else:
			for rxsplit in self.corso_fba.mapping[reaction]:
				corsofunction(rxsplit)

		mfunction(reaction)

	def find_reaction_limits(self, rx):
		self.corso_fba.cbmodel.set_objective({rx: 1}, True)
		fmin = self.corso_fba.cbmodel.optimize().objective_value()
		self.corso_fba.cbmodel.set_objective({rx: 1}, False)
		fmax = self.corso_fba.cbmodel.optimize().objective_value()

		return fmin, fmax

	def check_if_blocked(self, rx):
		fl = self.find_reaction_limits(rx)
		if fl[0] == nan and fl[1] == nan:
			return True
		else:
			return (abs(fl[0]) < 1e-6) and (abs(fl[1]) < 1e-6)

	def find_dependent_reactions(self, rx, constraint, constrainby, costfx, costbase, ntimes, eps):
		dependent, to_delete = self.__find_dependent_reactions(rx, constraint, constrainby, costfx, costbase, ntimes,
														   True,
														   eps)

		if self.lb[rx] < 0:
			bkw_dep, to_del_bkw = self.__find_dependent_reactions(rx, -constraint, constrainby, costfx, costbase,
																  ntimes,
																  False, eps)

			dependent = dependent | bkw_dep
			to_delete = to_del_bkw & to_delete

		return dependent, to_delete

	def __find_dependent_reactions(self, rx, constraint, constrainby, costfx, costbase, n_times, forward, eps):
		of_dict = {rx: 1}
		cost = costbase + costfx()
		# print(rx, cost)
		flux, corso_sol = self.corso_fba.optimize_corso(cost, of_dict, not forward, constraint, constrainby, eps=eps)

		dependent = abs(corso_sol.x()) > eps
		to_del = not dependent.any()
		if not to_del:
			for i in range(n_times - 1):
				cost = costbase + costfx()
				flux, corso_sol = self.corso_fba.optimize_corso(cost, of_dict, not forward, constraint, constrainby, eps=eps, flux1=flux)
				dependent = (abs(corso_sol.x()) > eps) | dependent
		else:
			dependent = zeros(dependent.shape).astype(bool)
		return dependent, to_del


