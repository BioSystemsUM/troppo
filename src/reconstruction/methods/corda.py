from numpy import array, zeros, sqrt, hstack, vstack, append, ones, where, floor, random, intersect1d

from cobamp.core.linear_systems import SteadyStateLinearSystem
from cobamp.core.optimization import LinearSystemOptimizer


class CORDA():
	def __init__(self, S, lb, ub, properties):
		self.S = array(S)
		self.lb, self.ub = array(lb), array(ub)
		self.properties = properties

	def step_one(self, S, lb, ub, es, pr, np, om, met_tests, current_reactions, constraint, constraintby, nl, ntimes):
		S, lb, ub = S.copy(), lb.copy(), ub.copy()
		m, n = self.S.shape
		costbase = zeros(n, )
		costbase[pr] = sqrt(om)
		costbase[np] = om

		es_fwd_del, es_bkw_del, pr_p, np_p = [zeros(dim,) for dim in [len(es), len(es), len(pr), len(np)]]
		hc_to_mc, hc_to_nc = zeros(len(es), len(pr)), zeros(len(es), len(np))

		met_test_ids = None

		# add test reactions to the S matrix as inactive irreversible reactions
		# when it's necessary to add them to the model, simply change the bounds
		# using met_test_ids as an

		if met_tests:
			met_test_ids = array([S.shape[1] + i for i in range(len(met_tests))])
			S = hstack([S, self.get_tests_stoichiometry(S.shape[0], met_tests)])
			append(lb, array([0 for _ in range(len(met_tests))]))
			append(ub, array([0 for _ in range(len(met_tests))]))
			append(costbase, array([0 for _ in range(len(met_tests))]))

		# TODO: missing the steps on lines 108-115

		# check for dependencies


		def update_dependency(x):
			active_idx = where(abs(x) > 1e-6)[0]
			for pres, ids in [(pr_p, pr), (np_p, np), (hc_to_mc, pr), (hc_to_nc, np)]:
				pres[intersect1d(active_idx, ids)] = 1

		def random_cost():
			return nl * floor(1e4 * random.rand(len(costbase), ))/10000

		# TODO: Shorten this!
		for i,rx in enumerate(es):
			x1, f = self.check_reaction_dependency(S, lb, ub, rx, random_cost(), constraint, constraintby, nl, minimize=False)
			if f > 1e-6:
				update_dependency(x1)
				for k in range(ntimes):
					cost_inc = costbase+random_cost()
					xn, f = self.check_reaction_dependency(S, lb, ub, rx, cost_inc, constraint, constraintby, nl, minimize=False)
					update_dependency(xn)
			else:
				es_fwd_del[i] = 1

			if lb[rx] < 0:
				xb, f = self.check_reaction_dependency(S, lb, ub, rx, costbase+random_cost(),
													   -constraint, constraintby, nl, minimize=True)
				if f > 1e-6:
					update_dependency(xb)
					for k in range(ntimes):
						cost_inc = costbase + random_cost()
						xbn, f = self.check_reaction_dependency(S, lb, ub, rx, cost_inc, constraint, constraintby, nl,
															   minimize=True)
						update_dependency(xbn)
				else:
					es_bkw_del[i] = 1

			else:
				es_bkw_del[i] = 1

		return _


	def check_reaction_dependency(self, S, lb, ub, rx, costfx, constraint, constraintby, nl, minimize):
		obj_fx = zeros(S.shape[1], )
		obj_fx[rx] = 1
		# TODO: Check what subl is - line 119
		return self.corsoFBA(S, lb, ub, of=obj_fx, constraint=constraint, constraintby=constraintby,
							 costas=costfx, nl=nl, minimize=minimize)

	def corsoFBA(self, S, lb, ub, of, costas, minimize=False, constraint=1, constraintby='val', eps=1e-6):
		constraint = abs(constraint) if constraintby == 'perc' else constraint
		zero_result = zeros(S.shape[1])
		lso, ss = self.get_linear_model(S, lb, ub, of, minimize)
		sol = lso.optimize()
		f1_f, f1_x = sol.objective_value(), sol.x()

		if abs(f1_f) < eps:
			return zero_result

		if constraintby == 'perc':
			f1_f = f1_x[of != 0] * (constraint / 100)
		elif constraintby == 'val':
			if (f1_f > constraint and minimize) or (f1_f < constraint and not minimize):
				return zero_result
			else:
				f1_f = constraint
		else:
			raise (ValueError('Invalid `constraintby` parameter. Must be either \'val\' or \'perc\''))

		if len(costas) == 1:
			costas = ones(S.shape[1],)
		if len(costas) == S.shape[1]:
			append(costas, costas)
		elif len(costas) != 2 * len(S.shape[1]):
			raise (ValueError('Invalid costs length'))

		n_orig = S.shape[1]
		rev_idx = where(lb < 0 & ub >= 0)[0]

		S_n, lb_n, ub_n, rx_mapping = self.make_irreversible(S, lb, ub)
		S_n = vstack([S_n, zeros(S_n.shape[0])])
		S_n[-1, :] = hstack([costas[:n_orig], costas[n_orig + rev_idx]])
		append(lb_n, array([1e20]))
		append(ub_n, array([0]))
		col = zeros(S_n.shape[0], )
		col[-1] = -1
		S_n = hstack([S_n, col])
		S_ov = zeros(S_n.shape[1])
		S_ov[-1] = 1

		obj_fluxes = where(of != 0)[0]
		for rx in obj_fluxes:
			lb_n[rx] = ub_n[rx] = f1_f[rx]
			if isinstance(rx_mapping[rx], tuple):
				lb_n[rx_mapping[rx][1]] = ub_n[rx_mapping[rx][1]] = 0

		lso2, ssls2 = self.get_linear_model(S_n, lb_n, ub_n, S_ov, True)
		flux2 = lso2.optimize()

		flux_x = flux2[:n_orig]
		flux_x[rev_idx] = flux_x[rev_idx] - flux2.x()[n_orig:-1]
		flux_x[abs(flux_x) < 1e-8] = 0

		return flux_x, f1_f

	def make_irreversible(self, S, lb, ub):
		## TODO: Order should be S, Srb
		irrev = array([i for i in range(self.S.shape[1]) if not (lb[i] < 0 and ub[i] > 0)])
		rev = array([i for i in range(self.S.shape[1]) if i not in irrev])
		Sr = S[:, rev]
		offset = S.shape[1]

		rx_mapping = {k: v if k in irrev else [v] for k, v in dict(zip(range(offset), range(offset))).items()}
		for i, rx in enumerate(rev):
			rx_mapping[rx].append(offset + i)
		rx_mapping = {k: tuple(v) if isinstance(v, list) else v for k, v in rx_mapping.items()}

		S_new = hstack([S, -Sr])
		nlb, nub = zeros(S_new.shape[1]), zeros(S_new.shape[1])
		for orig_rx, new_rx in rx_mapping.items():
			if isinstance(new_rx, tuple):
				nub[new_rx[0]] = abs(lb[orig_rx])
				nub[new_rx[1]] = ub[orig_rx]
			else:
				nlb[new_rx], nub[new_rx] = lb[orig_rx], ub[orig_rx]

		return S_new, nlb, nub, rx_mapping

	def get_linear_model(self, S, lb, ub, of=None, minimize=False):
		ss = SteadyStateLinearSystem(S, lb, ub, ['V' + str(i) for i in range(len(S.shape[1]))])
		lso = LinearSystemOptimizer(ss)
		if of:
			ss.set_objective(of, minimize)
		return lso, ss

	def get_tests_stoichiometry(self, row_shape, tests):
		S_test = zeros(row_shape, len(tests))
		S_test[(tests, array(range(len(tests))))] = 1
		return S_test
