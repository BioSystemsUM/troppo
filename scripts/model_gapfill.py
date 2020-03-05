import os
from json import JSONDecoder

import pandas as pd
from cobra.io import read_sbml_model

from cobra.flux_analysis.parsimonious import pfba
from cobra.core.reaction import Reaction
from troppo.tasks.core import TaskEvaluator, Task
from troppo.methods.gapfill.pathway_analysis import SubEFMGapfill

# paths to necessary files
# this is still hardcoded
ROOT_FOLDER = '../troppo/scripts/ccle_breast_cancer'

MODEL_PATH = os.path.join(ROOT_FOLDER, 'Recon3D_301.xml')  # SBML model path
DATA_PATH = os.path.join(ROOT_FOLDER, 'CCLE_breast_cancer_preprocessed.csv')  # tabular csv with expression
TASKS_PATH = os.path.join(ROOT_FOLDER, 'nl2019_tasks_r3d_compact.json')  # JSON containing metabolic tasks
CMODL_PATH = os.path.join(ROOT_FOLDER, 'Recon3D_301_consistent.xml')  # Consistent SBML model path - recommended!
TASK_RESULTS_PATH = os.path.join(ROOT_FOLDER, 'results',
                                 'r3d_compact_task_results_ccle_bc_new_nodrains_only_feas.json')  # task evaluation
CS_MODEL_DF_PATH = os.path.join(ROOT_FOLDER, 'results',
                                'r3d_compact_ccle_bc_fastcore.csv')  # context-specific models extracted from algos
MEDIA_DICT_PATH = os.path.join(ROOT_FOLDER, 'r3d_media_metabolites.json')

with open(MEDIA_DICT_PATH, 'r') as f:
	media_metabolites = JSONDecoder().decode(f.read())

fastcore_res_dict = pd.read_csv(CS_MODEL_DF_PATH, index_col=[0, 1]).T.to_dict()
# Task evaluation #

# read the original model to avoid compatibility issues with the tasks (e.g. missing metabolites from the block)
task_model = read_sbml_model(CMODL_PATH)
for k in task_model.boundary:
	k.knock_out()


# with task_model as mod:
# 	exchange = set(r.id for r in mod.boundary)
# 	all_rx = set(r.id for r in mod.reactions)
# 	missing = set([k for k,v in fastcore_res_dict[('fastcore', 'ACH-000019')].items() if not v and k not in exchange])-\
# 	          {'3DSPHR', 'ASPCTr', 'CBPS', 'RPI', 'SERPT', 'TAG_HSad_NE', 'TMDS', 'r0224'}
#
# 	for r in missing:
# 		mod.reactions.get_by_id(r).knock_out()
# 	for r in exchange:
# 		mod.reactions.get_by_id(r).bounds = (0, 1000)
#
#
#
# 	rpmi_dict = {k: 1000 for k in set(media_metabolites['DMEM']) | set(media_metabolites['RPMI'])
# 	             if k in [m.id for m in task_model.metabolites]}
# 	rpmi_dict.update({'o2[e]': 1000})
# 	for k in rpmi_dict:
# 		if 'EX_'+k in all_rx:
# 			mod.reactions.get_by_id('EX_'+k).bounds = (-1000, 1000)
# 		else:
# 			nrx = Reaction('EX_'+k, lower_bound=-1000, upper_bound=1000)
# 			nrx.add_metabolites({mod.metabolites.get_by_id(k): 1})
# 			mod.add_reactions([nrx])
# 	#mod.genes.get_by_id('10797.1').knock_out()
# 	mod.objective = 'biomass_reaction'
# 	sol = mod.optimize()
# 	min_sol = pfba(mod)
#
# print(task_model.summary(min_sol))
#

with open(TASK_RESULTS_PATH, 'r') as f:
	task_eval_results = JSONDecoder().decode(f.read())

# Gapfill model #
# This example will involve the creation of a task that evaluates whether a model can produce biomass on a given
# growth media. We will use one sample - ACH-000624 - grown on RPMI and FBS as medium
# Tasks associated with growth could be supplied, but in this example we will create them from scratch

# First, we load a previously generated JSON file with several commercially available medium formulations
# This file encodes a dictionary mapping media with the metabolites they add to the system

# We first need a dictionary with the inputs and their bounds. Since we can't assume every single compound is
# used, their bounds are left as [0, 1000]. Metabolites should also be present in the model
special = ['glc_D[e]']
rpmi_dict = {k: [0, 1000] if k not in special else [1,1000] for k in media_metabolites['MEM'] if k in [m.id for m in task_model.metabolites]}
rpmi_dict.update({'o2[e]': [0, 1000]}) # add oxygen uptake

# A partial task is then created, defining only the input metabolites - it is worth noting this task should actually
# fail on its own, but we will only use this as a component for the final task
rpmi_task = Task(should_fail=False, inflow_dict=rpmi_dict, name='RPMI_medium')

# The second component involves adding a biomass reaction with a forced bound - alternatively, a bound could be
# added on the existing biomass reaction, but this demonstrates the flexibility of this approach
# Several biomass compositions can be defined and added as tasks
# We first extract the generic biomass composition coefficients
biomass_eqn = {k.id: v for k, v in task_model.reactions.biomass_reaction.metabolites.items()}

# Since the medium should provide the necessary inputs for growth, we must define the outputs. This is a very simple
# definition that assumes no by-products. These can be added if needed
#biomass_outputs = {k: [0, 1000] for k, v in biomass_eqn.items() if v > 0}
biomass_outputs = {}
# To avoid infeasibility due to missing export reactions, we also add the remaining extracellular metabolites to the
# output dictionary, but with open bounds
biomass_outputs.update({k.id: [0, 1000] for k in task_model.metabolites if k.compartment == 'e' and
                        (k.id not in biomass_outputs.keys()) and (k.id not in rpmi_dict.keys())})

# We can now create the second sub-task, involving an added reaction representing growth with a minimum flux of 1
# and also defining the outputs of the biomass reaction
biomass_task = Task(should_fail=False, reaction_dict={'generic_biomass': (biomass_eqn, [1, 1000])}, name='growth',
                    outflow_dict=biomass_outputs)

# Both tasks can be combined with the addition operator, merging components from both
growth_task = biomass_task + rpmi_task
# The name is then changed to avoid whitespace characters
growth_task.name = 'growth_on_RPMI'

# We can now assess whether a generic human cell is capable of growing in the RPMI medium. This should not fail.
# Please note that boundary reactions should be removed before evaluating any tasks.
tev = TaskEvaluator(model=task_model, tasks=[growth_task], solver='CPLEX')
tev.current_task = growth_task.name
res = tev.evaluate()
print(res[0])

# This analysis can be extended to the results we also examined earlier
# rpmi_growth_results = {}
# for sample, context in fastcore_res_dict.items():
# 	def context_function(model):
# 		for rx in [k for k,v in context.items() if not v and k in model.reaction_names]:
# 			model.set_reaction_bounds(rx, lb=0, ub=0)
# 	rpmi_growth_results[sample] = tev.evaluate(context_function)[0]

# Upon inspecting the results, we can clearly see that none of the generated models is capable of growth.
# Gapfilling methods can be used to tackle this issue. Although any task can be corrected with gapfill, we will
# only try to gapfill until we can achieve growth, since it is the only solid assumption we can make about these
# cell lines.



task_reactions = set(range(len(task_model.reactions), len(tev.model.reaction_names))) | set([tev.model.reaction_names.index(k) for k in set([t.id for t in task_model.boundary])])
exchange = set([tev.model.reaction_names.index(k) for k in set([t.id for t in task_model.boundary])])
forced_reactions = set(
	[tev.model.decode_index(k, 'reaction') for k, v in growth_task.get_task_bounds().items() if '_inflow' in k])
forced_reactions -= set([k for k in forced_reactions if tev.model.optimize({k: 1}).objective_value() < 1e-10])

efm_gapfill = SubEFMGapfill(template_model=tev.model, task_reactions=task_reactions, solver='CPLEX',
                            subset_forced_reactions=set(), big_m=True, big_m_value=100000, max_threads=12)
gapfill_reaction_dict = {}
sample = ('fastcore', 'ACH-000019')

for sample in list(fastcore_res_dict.keys())[:1]:
	missing = {tev.model.map_labels['reaction'][k] for k, v in
	           fastcore_res_dict[sample].items() if (not v) and (k in tev.model.reaction_names)}

	gapfill_result = efm_gapfill.gapfill(missing_set=missing, forced=forced_reactions,
	                                     non_forced=task_reactions - forced_reactions - exchange)
	gapfill_reactions = {tev.model.reaction_names[f] for f in
	                     {k for k, v in gapfill_result[0].items() if abs(v) > 1e-8} & missing}
	gapfill_reaction_dict[sample] = gapfill_reactions


from cobra.core.reaction import Reaction

forced_bounds = {'gln_L[e]': 1.8709677,'arg_L[e]': 0.5971564,'cys_L[e]': 0.1,'his_L[e]': 0.2,'ile_L[e]': 0.39694658,
                 'leu_L[e]': 0.39694658, 'lys_L[e]': 0.3989071, 'met_L[e]': 0.10067114, 'phe_L[e]': 0.19393939,
                 'thr_L[e]': 0.40336135, 'trp_L[e]': 0.04901961, 'tyr_L[e]': 0.19889502, 'val_L[e]': 0.3931624,
                 'chol[e]': 0.007142857, 'pnto_R[e]': 0.002096436, 'fol[e]': 0.0022675737, 'ncam[e]': 0.008196721,
                 'pydx[e]': 0.004901961, 'ribflv[e]': 0.00026595744, 'thmtp[e]': 0.002967359, 'inost[e]': 0.011111111,
                 'cl[e]': 124.37063180000001, 'so4[e]': 0.8130081, 'k[e]': 5.3333335, 'hco3[e]': 26.190475,
                 'na1[e]': 118.25420050000001, 'h[e]': 1.0128205, 'pi[e]': 1.0128205,
                 'glc_D[e]': 5.5555553, 'ala_L[e]': 1.8709677, 'o2[e]': 1000
}

essentials = ['glc_D[e]']
all_rx = set(r.id for r in task_model.reactions)

manual_gpfl = set(['EX_h2o[e]'])
non_context = set([task_model.reactions.get_by_id(k) for k in (set(tev.model.reaction_names[m] for m in missing) -
                                                               set(gapfill_reactions | manual_gpfl))])
non_context_rxs = [r.id for r in non_context]
with task_model as mod:
	for r in non_context:
		r.knock_out()
	for k in mod.boundary:
		k.bounds = (0,1000)
	for k in set(rpmi_dict.keys()) | set(forced_bounds.keys()):
		if 'EX_'+k in all_rx:
			mod.reactions.get_by_id('EX_'+k).bounds = (-1, -0)
		else:
			nrx = Reaction('EX_'+k, lower_bound=-1, upper_bound=0)
			nrx.add_metabolites({mod.metabolites.get_by_id(k): 1})
			mod.add_reactions([nrx])
	for k,v in forced_bounds.items():
		mod.reactions.get_by_id('EX_'+k).bounds = (-1, 0)
	# for k in essentials:
	# 	mod.reactions.get_by_id('EX_'+k).upper_bound = -forced_bounds[k]/10
	mod.reactions.get_by_id('EX_h2o[e]').bounds = (-1000, 1000)
	# mod.reactions.get_by_id('EX_glc_D[e]').bounds = (-5.5555553, -5.5555553)
	# mod.reactions.get_by_id('EX_o2[e]').bounds = (-1000, -0.5)
	mod.objective = 'biomass_reaction'
	mod.reactions.biomass_reaction.bounds = (0, 1000)
	sol = pfba(mod)
	sol2 = mod.optimize()

print(task_model.summary(sol, threshold=0.000001))
print(task_model.summary(sol2, threshold=0.1))

sol.fluxes['PGI']
