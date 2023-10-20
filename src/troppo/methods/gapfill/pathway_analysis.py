from cobamp.core.linear_systems import IrreversibleLinearSystem
from cobamp.algorithms.kshortest import KShortestEnumerator
from cobamp.core.models import ConstraintBasedModel
from numpy import zeros, array, where

from collections import Counter


class SubEFMGapfill(object):
	"""
	Solve the gap-filling problem using an Elementary Flux Mode enumeration algorithm.

	Parameters
	----------
	template_model: dict
		The template model for the gap-filling problem. Must contain {S, lb, ub} and possibly rx_names/met_names.
	task_reactions
		Contains the task reactions.
	subset_forced_reactions
		Contains the indexes of the forced reactions' subset.
	iterative: bool
		Flag to determine whether the k-shortest algorithm behaviour.
	max_time: int
		Maximum time for the algorithm to run.
	max_threads: int
		Maximum number of threads to be used.
	big_m: bool
		Flag to determine whether to use big-M constraints.
	big_m_value: int
		Value of the big-M constraints.
	solver: str
		Solver to be used by the IrreversibleLinearPatternSystem.
	at_most_n_sols: int
		Maximum number of solutions to be returned.
	populate_max_size: int
		Maximum size of the population.
	"""
	def __init__(self, template_model: dict, task_reactions, subset_forced_reactions, iterative=True,
				 max_time: int = 0, max_threads: int = 0, big_m: bool = False, big_m_value: int = 1000,
				 solver: str = None, at_most_n_sols: int = 1, populate_max_size: int = None):
		# media: dict - must contain the following keys - consumed, non_consumed, produced.
		# Values can be set as empty lists if not needed

		if isinstance(template_model, dict):
			S, lb, ub = [template_model[k] for k in ['S', 'lb', 'ub']]
			self.cb_model = ConstraintBasedModel(S=S, thermodynamic_constraints=list(zip(lb, ub)))
		elif isinstance(template_model, ConstraintBasedModel):
			self.cb_model = template_model

		self.cb_model.make_irreversible()
		self.task_reactions = task_reactions
		self.subset_forced_reactions = subset_forced_reactions
		self.__iterative = iterative
		self.__max_time = max_time
		self.__max_threads = max_threads
		self.__big_m = big_m
		self.__big_m_value = big_m_value
		self.__at_most_n_sols = at_most_n_sols
		self.__popmaxsz = populate_max_size
		self.enumerated_solutions = []

		Si = self.cb_model.get_stoichiometric_matrix()
		lbi, ubi = list(zip(*self.cb_model.bounds))

		lsys = IrreversibleLinearSystem(
			S=Si, lb=array(lbi), ub=array(ubi), solver=solver)

		self.enumerator_obj = KShortestEnumerator(
			linear_system=lsys,
			m_value=self.__big_m_value,
			force_big_m=self.__big_m,
			n_threads=self.__max_threads,
			max_time=self.__max_time,
		)

		self.__mask_len = Si.shape[1]

	def gapfill(self, missing_set, forced=(), non_forced=()):
		"""
		Gapfill the model using the task reactions and the subset of forced reactions.

		Parameters
		----------
		missing_set
			Dictionary containing the indexes of the missing reactions.
		forced:
			Dictionary containing the indexes of the forced reactions.
		non_forced:
			Dictionary containing the indexes of the non-forced reactions.

		Returns
		-------
		List of solutions.
		"""

		self.enumerator_obj.reset_enumerator_state()
		dvmap = self.enumerator_obj.get_model().get_dvar_mapping()

		final_subset = missing_set | set(self.subset_forced_reactions)
		subset_mask = zeros(self.__mask_len)
		subset_mask[list(final_subset)] = 1

		# forced_on = task reactions
		# forced_off = non task reactions
		n_ivars = len(self.enumerator_obj.model.get_dvars())

		non_forced_rxs = set(non_forced)
		non_forced_rxs |= set([dvmap[i][1] for i in non_forced if isinstance(dvmap[i], (list, tuple))])

		io_rxs = set(forced)
		io_rxs |= set([dvmap[i][1] for i in io_rxs if isinstance(dvmap[i], (list, tuple))])
		non_io_rxs = set(self.task_reactions)
		non_io_rxs |= set([dvmap[i][1] for i in non_io_rxs if isinstance(dvmap[i], (list, tuple))])

		non_io_rxs -= io_rxs
		io_rxs_mask, non_io_rxs_mask = zeros(n_ivars, ), zeros(n_ivars, )

		io_mask_idx = array(list([int(i) for i in (io_rxs - non_forced_rxs)]))
		non_io_mask_idx = array(list([int(i) for i in (non_io_rxs - non_forced_rxs)]))

		if len(io_mask_idx):
			io_rxs_mask[io_mask_idx] = 1

		if len(non_io_mask_idx):
			non_io_rxs_mask[non_io_mask_idx] = 1

		self.enumerator_obj.set_indicator_activity(forced_on=io_rxs_mask, forced_off=non_io_rxs_mask)
		self.enumerator_obj.set_objective_expression(subset_mask)

		if self.__iterative:
			enumerator = self.enumerator_obj.solution_iterator()
		else:
			enumerator = self.enumerator_obj.population_iterator(self.__popmaxsz)

		solutions = []
		has_next = True
		while has_next and (len(solutions) < self.__at_most_n_sols):
			try:
				sol_result = next(enumerator)
				if isinstance(sol_result, (list, tuple)):
					solutions.extend(sol_result)
				else:
					solutions.append(sol_result)
			except StopIteration as e:
				has_next = False
		self.enumerated_solutions = solutions
		return [sol.attribute_value(sol.SIGNED_INDICATOR_SUM) for sol in solutions]


class CombinatorialEFMGapfill():
	def __init__(self, template_model, media):  # tasks, media, min_acceptable_tasks=0.5):

		self.template = template_model
		# self.pass_fx = pass_fx
		if isinstance(self.template, ConstraintBasedModel):
			self.S = self.template.get_stoichiometric_matrix()
			self.lb, self.ub = map(array, self.template.get_bounds_as_list())
		elif isinstance(self.template, dict):
			self.S, self.lb, self.ub = [array(self.template[k]) for k in ['S', 'lb', 'ub']]

		self.cbmodel = ConstraintBasedModel(self.S, list(zip(self.lb, self.ub)))
		# self.tasks = tasks
		self.media = media

	# if isinstance(min_acceptable_tasks, int):
	# 	self.__min_tasks = min_acceptable_tasks
	# elif isinstance(min_acceptable_tasks, float) and (0 < min_acceptable_tasks <= 1):
	# 	self.__min_tasks = (min_acceptable_tasks // len(self.tasks)) if (len(self.tasks) > 0) else 0
	# else:
	# 	raise ValueError('Invalid parameter `min_acceptable_tasks`')
	def combinatorial_gapfill(self, sample_models, model_approval_fx):
		partials = self.generate_partial_models(sample_models)
		i = len(partials)
		gfr, success = {}, False
		while i >= 2 and not success:
			m1, m0 = partials[i - 2:i]
			gfr, success = self.gapfill_partial_pair(m0, m1, model_approval_fx)
			i -= 1
		return gfr, success

	def generate_partial_models(self, sample_models):
		c = Counter()
		for mod in sample_models:
			c.update(mod)

		partials = [set([k for k, v in c.items() if v >= n]) for n in range(len(sample_models) + 1)]
		return partials

	def prune_model(self, to_keep):
		tkid = list(to_keep)
		S, lb, ub = self.S.copy()[:, tkid], self.lb.copy()[tkid], self.ub.copy()[tkid]
		metabs = where(~(S == 0).all(axis=1))[0]
		S = S[metabs, :]
		return S, lb, ub, metabs

	def gapfill_partial_pair(self, m0, m1, model_approval_fx):
		assert len(m0) <= len(m1), 'Model 1 has less reactions than Model 0!'

		Sp, lbp, ubp, metabs = self.prune_model(m1)
		mapback = dict(zip(range(len(m1)), list(m1)))
		mapback_rev = dict(zip(list(m1), range(len(m1))))
		metab_mapback = dict(zip(list(metabs), range(len(metabs))))

		valid_problem = False
		try:
			pruned_media = {k: [metab_mapback[i] for i in v] for k, v in self.media.items()}
			valid_problem = True
		except:
			print('Partial model pair is unsuitable for the selected medium')

		if valid_problem:
			print('Attempting to gapfill partial model with', len(m0), 'reactions using', len(m1), 'sized template.')
			algorithm = SubEFMGapfill({'S': Sp, 'lb': lbp, 'ub': ubp}, subset_forced_reactions={},
									  media=pruned_media, big_m=True, big_m_value=1000)
			gfs = algorithm.gapfill(missing_set={mapback_rev[i] for i in list(m1 - m0)})
			if len(gfs) > 0:
				print('Found a solution!')
				gfrx = set([k for k, v in gfs[0].items() if abs(v) > 1e-9])
				tentative_m = m0 | {mapback[g] for g in gfrx}
				print('Original model has', len(m0), 'while the gapfilled solution has', len(tentative_m))
				is_valid_model = model_approval_fx(tentative_m)
				if is_valid_model:
					return tentative_m, True
				else:
					return {}, False
			else:
				return {}, False
		else:
			return {}, False


def simulate_context(model, context, objective, sense, sol_approval_fx):
	ko_count = 0
	ndim = range(model.get_stoichiometric_matrix().shape[1])
	for r in ndim:
		if r not in context:
			model.set_reaction_bounds(r, lb=0, ub=0, temporary=True)
			ko_count += 1
	try:
		print('Sanity check:', ko_count, 'knockouts applied.', len(set(ndim) - context), 'found.')
		real_ko_count = 0
		for i in ndim:
			lb, ub = model.get_reaction_bounds(i)
			if abs(lb) < 1e-9 and abs(ub) < 1e-9:
				real_ko_count += 1
		print('Sanity check:', real_ko_count, 'knockouts found within the CB Model')
		sol = model.optimize(coef_dict=objective, minimize=sense)
		print('Simulating model with', len(context), 'reactions...', sol.objective_value(), sol.status(),
			  sol_approval_fx(sol))
		return sol_approval_fx(sol)
	except:
		pass
	finally:
		model.revert_to_original_bounds()


if __name__ == '__main__':
	from urllib.request import urlretrieve
	from cobra.io import read_sbml_model
	from cobra.core import Reaction, Metabolite
	from cobamp.wrappers.cobra import COBRAModelObjectReader
	from random import sample, randint

	BIOMASS_RX_NAME = 'BIOMASS_Ecoli_core_w_GAM'
	BIOMASS_TRANS_NAME = 'BIOMASSt'
	BIOMASS_DRAIN_NAME = 'EX_biomass_e'
	BIOMASS_CYT_NAME = 'biomass_c'
	BIOMASS_EXT_NAME = 'biomass_e'


	def get_ecoli_core_model():

		model_url = 'http://bigg.ucsd.edu/static/models/e_coli_core.xml'
		model_path, _ = urlretrieve(model_url)
		ecoli_model = read_sbml_model(model_path)
		b_c, b_e = (Metabolite(BIOMASS_CYT_NAME, name='biomass (cytosol)'),
					Metabolite(BIOMASS_EXT_NAME, name='biomass (extracellular)'))

		ecoli_model.reactions.get_by_id(BIOMASS_RX_NAME).add_metabolites({b_c: 1})

		b_trans = Reaction(BIOMASS_TRANS_NAME, name='Biomass transport')
		b_trans.add_metabolites({b_c: -1, b_e: 1})

		b_drain = Reaction(BIOMASS_DRAIN_NAME, name='Biomass drain')
		b_drain.add_metabolites({b_e: -1})

		ecoli_model.add_reactions([b_trans, b_drain])

		return ecoli_model


	def random_knockouts(model_template, biomass_reactions_set, drains_set):
		exclude = biomass_reactions_set | drains_set
		choices = set(r.id for r in model_template.reactions) - set(exclude)
		missing_reaction_set = set(sample(choices, k=randint(len(choices) // 4, (3 * len(choices)) // 4)))
		return missing_reaction_set


	def simulate_model_with_kos(model_template, kos):
		with model_template as m:
			for rx in kos:
				m.reactions.get_by_id(rx).knock_out()
			sol = m.optimize()
		return sol


	model = get_ecoli_core_model()

	objreader = COBRAModelObjectReader(model)

	# s_matrix = objreader.get_stoichiometric_matrix()
	# lower_bounds, upper_bounds = objreader.get_model_bounds(separate_list=True)
	# template_model = {'S': s_matrix, 'lb': lower_bounds, 'ub': upper_bounds}
	template_model = objreader.to_cobamp_cbm('CPLEX')
	media = {'produced': [template_model.metabolite_names.index(BIOMASS_EXT_NAME)], 'non_consumed': [], 'consumed': []}

	drains = set([r.id for r in model.boundary])
	biomass_reactions = {BIOMASS_RX_NAME, BIOMASS_TRANS_NAME, BIOMASS_DRAIN_NAME}

	missing_set = random_knockouts(model, biomass_reactions, drains)
	final_subset = missing_set | biomass_reactions

	before_gapfill = model.optimize()

	missing_reactions = random_knockouts(model, biomass_reactions, drains)
	print(len(missing_reactions), 'total reactions missing.')

	after_missing_reactions = simulate_model_with_kos(model, missing_reactions)
	print('After removing selected KOs:', after_missing_reactions.status, after_missing_reactions.objective_value)

	missing_reaction_ids = [i for i, r_id in enumerate(template_model.reaction_names) if r_id in missing_set]
	biomass_reaction_ids = [template_model.reaction_names.index(k) for k in biomass_reactions]
	drain_reaction_ids = [template_model.reaction_names.index(k) for k in drains]

	gapfiller = SubEFMGapfill(template_model, biomass_reaction_ids, media, iterative=True, big_m=True, big_m_value=1000)
	gapfill_reactions = gapfiller.gapfill(missing_set=set(missing_reaction_ids), media=set(drain_reaction_ids))
	gapfill_rx_ids = set([template_model.reaction_names[i] for i, v in gapfill_reactions[0].items() if abs(v) > 1e-6])
	added_rx = set(missing_reactions) & gapfill_rx_ids
	print(len(added_rx), 'gapfilled reactions:', ','.join(added_rx))

	after_gapfilling = simulate_model_with_kos(model, set(missing_reactions) - gapfill_rx_ids)
	print('After gapfilling with EFM:', after_gapfilling.status, after_gapfilling.objective_value)

# any_sols = []
# sub_ko = None
