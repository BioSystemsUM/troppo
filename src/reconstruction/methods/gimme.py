import numpy as np

from cobamp.core.linear_systems import SteadyStateLinearSystem
from cobamp.core.optimization import LinearSystemOptimizer


class GIMME():
	def __init__(self, S, lb, ub, properties):
		self.S = np.array(S)
		self.lb, self.ub = np.array(lb), np.array(ub)
		self.properties = properties

	def preprocess(self, S, lb, ub, exp_vector, objectives):
		## make irreversible model
		Sn, lb_n, ub_n, irrev_mapping = self.make_irreversible(S, lb, ub)

		## check for expression data size
		exp_vector_n = np.zeros(Sn.shape[1], )
		for rxn, val in enumerate(exp_vector):
			rxmap = irrev_mapping[rxn]
			if isinstance(rxmap, tuple):
				exp_vector_n[rxmap[0]] = exp_vector_n[rxmap[1]] = val
			else:
				exp_vector_n[rxmap] = val

		## adjust objectives for irreversible model
		objectives = [self.adjust_objective_to_irreversible(Sn, obj, irrev_mapping) for obj in objectives]

		return Sn, np.array(lb_n), np.array(ub_n), irrev_mapping, exp_vector_n, objectives

	def run(self):

		S, lb, ub, exp_vector, objectives = self.S, self.lb, self.ub, self.properties['exp_vector'], self.properties[
			'objectives']

		if self.properties['preprocess']:
			S, lb, ub, irrev_mapping, exp_vector, objectives = self.preprocess(S, lb, ub, exp_vector, objectives)

		gimme_solution = self.solve_gimme_problem(S, lb, ub, exp_vector, objectives)

		if self.properties['preprocess']:
			gimme_solution = np.array(
				[np.max(gimme_solution[np.array(new)]) if isinstance(new, tuple) else gimme_solution[new] for orig, new
				 in irrev_mapping.items()])

		return gimme_solution

	def solve_gimme_problem(self, S, lb, ub, exp_vector, objectives):

		obj_frac = self.properties['obj_frac']
		flux_thres = self.properties['flux_threshold']

		M, N = S.shape
		## FBA for each objective
		lsystem = SteadyStateLinearSystem(S, lb, ub, ['V' + str(i) for i in range(S.shape[1])])
		lso = LinearSystemOptimizer(lsystem)

		def find_objective_value(obj):
			lsystem.set_objective(obj, False)
			sol = lso.optimize()
			return sol.objective_value()

		objective_values = list(map(find_objective_value, objectives))

		gimme_model_objective = np.array(
			[flux_thres - exp_vector[i] if -1 < exp_vector[i] < flux_thres else 0 for i in range(N)])

		objective_lbs = sum(list(map(lambda a, b: a * b * obj_frac, objectives, objective_values)))
		objective_ids = np.nonzero(objective_lbs)[0]
		lb[objective_ids] = objective_lbs[objective_ids]

		gimme_system = SteadyStateLinearSystem(S, lb, ub, var_names=['GIMME' + str(i) for i in range(N)])
		gimme_lso = LinearSystemOptimizer(gimme_system)
		gimme_system.set_objective(gimme_model_objective, True)

		gimme_solution = gimme_lso.optimize()

		return self.get_reaction_activity(gimme_solution, exp_vector, flux_thres)

	def adjust_objective_to_irreversible(self, S_new, objective, mapping):
		objective_new = np.zeros(S_new.shape[1], )
		nzids = np.nonzero(objective)[0]
		for id in nzids:
			objective_new[np.array(mapping[id])] = objective[id]
		return objective_new

	def make_irreversible(self, S, lb, ub):
		## TODO: Order should be S, Srb
		irrev = np.array([i for i in range(self.S.shape[1]) if not (lb[i] < 0 and ub[i] > 0)])
		rev = np.array([i for i in range(self.S.shape[1]) if i not in irrev])
		Si, Sr = S[:, irrev], S[:, rev]
		offset = Si.shape[1]
		rx_mapping = dict(zip(irrev, range(offset)))
		rx_mapping.update({ii: (offset + n, offset + n + Sr.shape[1]) for n, ii in
						   enumerate([i for i in range(S.shape[1]) if i not in irrev])})

		S_new = np.hstack([Si, Sr, -Sr])
		nlb, nub = np.zeros(S_new.shape[1]), np.zeros(S_new.shape[1])
		for orig_rx, new_rx in rx_mapping.items():
			if isinstance(new_rx, tuple):
				nub[new_rx[0]] = abs(lb[orig_rx])
				nub[new_rx[1]] = ub[orig_rx]
			else:
				nlb[new_rx], nub[new_rx] = lb[orig_rx], ub[orig_rx]

		return S_new, nlb, nub, rx_mapping

	def get_reaction_activity(self, solution, exp_vector, flux_threshold):
		gimme_fluxes = np.array([kv[1] for i, kv in enumerate(solution.var_values().items())])
		activity = np.zeros(gimme_fluxes.shape)
		ones = (exp_vector > flux_threshold) | (exp_vector == -1)
		twos = gimme_fluxes > 0
		activity[ones] = 1
		activity[twos & ~ones] = 2

		return activity
