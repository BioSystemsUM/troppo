from cobamp.wrappers.method_wrappers import KShortestEFMEnumeratorWrapper
from cobamp.core.models import ConstraintBasedModel

class SubEFMGapfill(object):
	def __init__(self, template_model, subset_forced_reactions, media, iterative=True):
		'''
		:param template_model: dict - must contain {S, lb, ub} and possibly rx_names/met_names
		:param subset_forced_reactions: list - contains reaction indexes
		:param media: dict - must contain the following keys - consumed, non_consumed, produced
		values can be set as empty lists if not needed
		:param iterative: boolean flag - determines whether the k-shortest algorithm behaviour
		'''
		if isinstance(template_model, dict):
			S, lb, ub = [template_model[k] for k in ['S', 'lb', 'ub']]
			self.cb_model = ConstraintBasedModel(S=self.S, thermodynamic_constraints=list(zip(lb, ub)))
		elif isinstance(template_model, ConstraintBasedModel):
			self.cb_model = template_model

		self.subset_forced_reactions = subset_forced_reactions
		self.met_pr, self.met_co, self.met_nc =  [media[k] for k in ['produced','consumed','non_consumed']]
		self.__iterative = iterative

	def gapfill(self, missing_set, at_most_n_sols=1):
		algo_type = 'kse_iterative' if self.__iterative else 'kse_populate'
		final_subset = missing_set | self.subset_forced_reactions
		stop_criteria = at_most_n_sols if self.__iterative else len(final_subset)
		algorithm = KShortestEFMEnumeratorWrapper(self.cb_model, subset=final_subset, stop_criteria=stop_criteria,
			non_consumed=self.met_nc, consumed=self.met_co, algorithm_type=algo_type, produced=self.met_pr)

		enumerator = algorithm.get_enumerator()
		solutions = []
		while len(solutions) < at_most_n_sols:
			try:
				sols = next(enumerator)
				if self.__iterative:
					solutions.append(sols)
				else:
					solutions.extend(sols)
			except StopIteration as e:
				print('No more solutions')

		return [set(d.keys()) for d in solutions]

if __name__ == '__main__':
	from urllib.request import urlretrieve
	from cobra.io.sbml3 import read_sbml_model
	from cobra.core import Reaction, Metabolite
	from cobamp.wrappers import COBRAModelObjectReader, KShortestEFMEnumeratorWrapper
	from random import sample, randint

	BIOMASS_RX_NAME = 'BIOMASS_Ecoli_core_w_GAM'
	BIOMASS_TRANS_NAME = 'BIOMASSt'
	BIOMASS_DRAIN_NAME = 'EX_biomass_e'
	BIOMASS_CYT_NAME = 'biomass_c'
	BIOMASS_EXT_NAME = 'biomass_e'

	def get_ecoli_core_model():


		model_url = 'http://bigg.ucsd.edu/static/models/e_coli_core.xml'
		model_path, _ = urlretrieve(model_url)
		model = read_sbml_model(model_url)
		b_c, b_e = (Metabolite(BIOMASS_CYT_NAME, name='biomass (cytosol)'),
					Metabolite(BIOMASS_EXT_NAME, name='biomass (extracellular)'))

		model.reactions.get_by_id(BIOMASS_RX_NAME).add_metabolites({b_c: 1})

		b_trans = Reaction(BIOMASS_TRANS_NAME, name='Biomass transport')
		b_trans.add_metabolites({b_c: -1, b_e: 1})

		b_drain = Reaction(BIOMASS_DRAIN_NAME, name='Biomass drain')
		b_drain.add_metabolites({b_e: -1})

		model.add_reaction(b_trans)
		model.add_reaction(b_drain)

		return model


	def random_knockouts(model, biomass_reactions, drains):
		exclude = biomass_reactions | drains
		choices = set(r.id for r in model.reactions) - set(exclude)
		missing_set = set(sample(choices, k=randint(len(choices) // 4, (3 * len(choices)) // 4)))
		return missing_set

	def simulate_model_with_kos(model, kos):
		with model as m:
			for rx in kos:
				m.reactions.get_by_id(rx).knock_out()
			sol = m.optimize()
		return sol


	model = get_ecoli_core_model()

	objreader = COBRAModelObjectReader(model)

	S = objreader.get_stoichiometric_matrix()
	lb, ub = objreader.get_model_bounds(separate_list=True)
	#template_model = {'S': S, 'lb': lb, 'ub': ub}
	template_model = objreader.to_cobamp_cbm('CPLEX')
	media = {'produced':[BIOMASS_EXT_NAME], 'non_consumed':[], 'consumed':[]}

	drains = set([r.id for r in model.boundary])
	biomass_reactions = {BIOMASS_RX_NAME, BIOMASS_TRANS_NAME, BIOMASS_DRAIN_NAME}

	missing_set = random_knockouts(model, biomass_reactions, drains)
	final_subset = missing_set | biomass_reactions

	before_gapfill = model.optimize()

	missing_reactions = random_knockouts(model, biomass_reactions, drains)
	print(len(missing_reactions),'total reactions missing.')

	after_missing_reactions = simulate_model_with_kos(model, missing_reactions)
	print('After removing selected KOs:',after_missing_reactions.status,after_missing_reactions.objective_value)

	gapfiller = SubEFMGapfill(template_model, biomass_reactions, media)
	gapfill_reactions = gapfiller.gapfill(missing_set=set(missing_reactions))
	print(len(gapfill_reactions[0]),'gapfilled reactions:',','.join(gapfill_reactions[0]))

	after_gapfilling = simulate_model_with_kos(model, set(missing_reactions) - gapfill_reactions[0])
	print('After gapfilling with EFM:',after_gapfilling.status,after_gapfilling.objective_value)

# any_sols = []
# sub_ko = None

