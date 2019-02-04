from numpy import array, zeros, sqrt, vstack, where, floor, random, unique, \
	apply_along_axis, nan, empty

from cobamp.core.models import CORSOModel
from cobamp.core.models import ConstraintBasedModel


class CORDA():
	def costfx_factory(self, nl, om, costbase):
		return lambda: nl * floor(om * random.rand(len(costbase), )) / om

	def __init__(self, S, lb, ub, properties):
		self.S = array(S)
		self.lb, self.ub = array(lb), array(ub)
		self.properties = properties

		rx_names, mt_names = ['V' + str(i) for i in range(S.shape[1])], ['M' + str(i) for i in range(S.shape[0])]
		cbmodel = ConstraintBasedModel(S, list(zip(lb, ub)), reaction_names=rx_names,
									   metabolite_names=mt_names, optimizer=True)
		self._m, self._n = self.S.shape
		self.corso_fba = CORSOModel(cbmodel)

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

		return self.step_two(rx_cat, constraint, constrainby, nl, ntimes, om, pr_to_np)


	def step_two(self, rx_cat, constraint, constrainby, nl, ntimes, om=1e4, pr_to_np=2):
		rx_cat = array(rx_cat)
		costbase = zeros(self._n, )

		costbase[rx_cat == 2] = sqrt(om)
		costbase[rx_cat == 3] = om

		HC_reactions = where(rx_cat == 1)[0]
		costfx = self.costfx_factory(nl, om, costbase)
		for rx in HC_reactions:
			dep, to_del = self.find_dependent_reactions(rx, constraint, constrainby, costfx, costbase, ntimes, eps=1e-6)
			rx_cat[dep] = 1
			if to_del:
				rx_cat[rx] = -1

		costbase = zeros(self._n, )
		costbase[rx_cat == 3] = om

		PR_reactions = where(rx_cat == 2)[0]
		NP_reactions = where(rx_cat == 3)[0]

		costfx = self.costfx_factory(nl, om, costbase)

		PR_NP = []
		for i, rx in enumerate(PR_reactions):
			dep, to_del = self.find_dependent_reactions(rx, constraint, constrainby, costfx, costbase, ntimes, eps=1e-6)
			PR_NP.append(dep[NP_reactions])
			if to_del:
				rx_cat[rx] = -1

		PR_NP = vstack(PR_NP)
		NP_occurrence = apply_along_axis(sum, 0, PR_NP)
		np_to_pr_idx = NP_reactions[NP_occurrence > pr_to_np]

		if len(np_to_pr_idx) > 0:
			rx_cat[np_to_pr_idx] = 2
			PR_NP = PR_NP[:, NP_occurrence > pr_to_np]
			if PR_NP.shape[1] > 0:
				PR_NP = vstack([zeros(len(np_to_pr_idx), PR_NP.shape[1])])
			else:
				PR_NP = array([])

		# 2.2
		PR_reactions = where(rx_cat == 2)[0]
		NP_reactions = where(rx_cat == 3)[0]

		PR_to_check_l8r = []
		for rx in NP_reactions:
			if isinstance(self.corso_fba.mapping[rx], int):
				self.corso_fba.set_reaction_bounds(rx, lb=0, ub=0)
			else:
				for rxsplit in self.corso_fba.mapping[rx]:
					self.corso_fba.set_reaction_bounds(rxsplit, lb=0, ub=0)

			self.corso_fba.cbmodel.set_reaction_bounds(rx, lb=0, ub=0)

		res2 = []
		for i, rx in enumerate(PR_reactions):
			to_del = self.check_if_blocked(rx)
			if to_del:
				PR_to_check_l8r.append(rx)
				if len(np_to_pr_idx) > 0:
					np_from_rx = where(PR_NP[i, :] > 0)[0]
					if len(np_from_rx) == 0:
						print('Undefined')
					else:
						for kn in np_from_rx:
							res2.append(kn)
				rx_cat[rx] = -1
		res2 = unique(sorted(array(res2)))

		PR_reactions = where(rx_cat == 2)[0]
		rx_cat[PR_reactions] = 1
		rescued = set(sorted(PR_to_check_l8r + res2.tolist()))

		ES_reactions = rx_cat == 1
		OT_reactions = rx_cat == 0

		to_block = where(ES_reactions | OT_reactions)[0]

		block_corso = lambda rx: self.corso_fba.set_reaction_bounds(rx, lb=0, ub=0)
		block_cbmodel = lambda rx: self.corso_fba.cbmodel.set_reaction_bounds(rx, lb=0, ub=0)

		for rx in to_block:
			self.do_function_for_reactions_on_both_models(rx, block_cbmodel, block_corso)

		costbase = zeros(self._n, )
		costbase[OT_reactions] = om

		ES_OT = []
		for rx in where(ES_reactions)[0]:
			dep, to_del = self.find_dependent_reactions(rx, constraint, constrainby, costfx, costbase, ntimes, 1e-6)
			ES_OT.append(dep[OT_reactions])

		OT_reaction_ids = where(OT_reactions)[0]
		ES_OT = vstack(ES_OT)
		if ES_OT.shape[1] > 0:
			OT_occurrence = apply_along_axis(sum, 0, ES_OT)
			ot_to_es_idx = OT_reaction_ids[OT_occurrence != 0]
			rx_cat[ot_to_es_idx] = 1

		return rx_cat

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

			dependent = dependent & bkw_dep
			to_delete = to_del_bkw & to_delete

		return dependent, to_delete

	def __find_dependent_reactions(self, rx, constraint, constrainby, costfx, costbase, n_times, forward, eps):
		of_dict = {rx: 1}
		cost = costbase + costfx()
		_, corso_sol = self.corso_fba.optimize_corso(cost, of_dict, not forward, constraint, constrainby)

		dependent = corso_sol.x() > eps
		to_del = not dependent.any()
		if not to_del:
			for i in range(n_times - 1):
				cost = costbase + costfx()
				_, corso_sol = self.corso_fba.optimize_corso(cost, of_dict, not forward, constraint, constrainby)
				dependent = (corso_sol.x() > eps) | dependent
		else:
			dependent = zeros(dependent.shape).astype(bool)
		return dependent, to_del
