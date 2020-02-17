from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.tasks.core import TaskEvaluator
from troppo.tasks.task_io import JSONTaskIO
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from cobra.flux_analysis.variability import find_blocked_reactions
from cobamp.utilities.parallel import batch_run
import pandas as pd
from json import JSONEncoder

from numpy import log
import re

if __name__ == '__main__':

	# helper functions
	patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')  # find .{number} references
	replace_alt_transcripts = lambda x: patt.sub('', x)  # replace .{number} with nothing

	import os

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

	# Context-specific model reconstruction #

	# model preprocessing only if the model isn't loaded already
	# this consists of:
	# - removing artificial sinks and drug modules
	# - removing blocked reactions
	if not os.path.exists(CMODL_PATH):
		model_consistent = read_sbml_model(MODEL_PATH)
		model_consistent.remove_reactions(
			[r for r in model_consistent.reactions if r.id[:3] == 'DM_' or r.id[:5] == 'sink_'], remove_orphans=True)
		blocked_reactions = find_blocked_reactions(model_consistent)
		model_consistent.remove_reactions(blocked_reactions, remove_orphans=True)
		write_sbml_model(model_consistent, CMODL_PATH)  # write a model file if it doesn't exist
	else:
		model_consistent = read_sbml_model(CMODL_PATH)

	# read the csv on DATA_PATH as a transcriptomics data set and get OmicsContainers, one per sample
	ocs = TabularReader(path=DATA_PATH, nomenclature='entrez_id', omics_type='transcriptomics').to_containers()

	# create a reconstruction wrapper instance that loads the consistent model
	# GPRs are assumed to be in DNF form if the ratio between tokens and unique genes is above ttg_ratio
	# the helper function replace_alt_transcripts is used here to replace .{number} with nothing
	rw = ReconstructionWrapper(model_consistent, ttg_ratio=9999, gpr_gene_parse_function=replace_alt_transcripts)

	# gene score threshold to determine core reactions
	# on this dataset, genes were normalized using the following rule : 5*log(1+(expression/mean))
	# thus, if the expression is average, the final score will be 5*log(2)
	# we intend to keep reactions above this threshold as core
	t = (5 * log(2))


	# since we will be running multiple samples, we can generalize a model reconstruction as a function
	# this function should take two arguments
	# the first will be the variable between samples (in this case, a different omics container)
	# the second should be a dictionary with static variables needed to build the model:
	# - threshold
	# - reconstruction wrapper
	def reconstruction_func(omics_container, params):
		t, rw = [params[k] for k in ['t', 'rw']]  # load parameters
		try:
			# if no errors appear, call the run_from_omics method passing the omics_container,
			# algorithm string, integration strategy (how the core is determined) and a solver
			# for fastcore, a threshold-based integration strategy retrieves core reactions if the score
			# is above the threshold t
			return rw.run_from_omics(omics_container=omics_container, algorithm='fastcore',
									 integration_strategy=('threshold', t), solver='CPLEX')
		except:
			# the result from run_from_omics is a dict mapping reaction ids and a boolean flag - True if
			# the reaction is in the model or false otherwise
			# in case an error arises, assume all reactions are False
			return {r: False for r in rw.model_reader.r_ids}


	# parallel reconstruction can be achieved with the batch_run function that takes in 4 key arguments:
	# - the function to apply multiple times - reconstruction_func
	# - a list of objects for which we want to apply the function - ocs (list of omics_containers)
	# - a dictionary with the static params - containing a 't' and 'rw' entry (see reconstruction_func)
	# - an integer value specifying the amount of parallel processes to be run
	batch_fastcore_res = batch_run(reconstruction_func, ocs, {'t': t, 'rw': rw}, threads=11)

	# create a dict mapping tuple of algorithm and condition informations: e.g. ('fastcore','MCF7')
	# to each result from fastcore
	fastcore_res_dict = dict(zip([('fastcore', oc.condition) for oc in ocs], batch_fastcore_res))

	# write these results as a dataframe for future reference
	pd.DataFrame.from_dict(fastcore_res_dict, orient='index').to_csv(CS_MODEL_DF_PATH)

	# Task evaluation #

	# parse tasks from a previously existing JSON
	# the supplied file contains tasks adapted from the publication of Richelle et. al, 2019
	task_list = JSONTaskIO().read_task(TASKS_PATH)

	# read the original model to avoid compatibility issues with the tasks (e.g. missing metabolites from the block)
	task_model = read_sbml_model(MODEL_PATH)

	# tasks should be evaluated without open boundary reactions. we can easily close them on the COBRA model
	for k in task_model.boundary:
		k.knock_out()

	# get the names of all reactions in the model - this will be useful further on
	all_reactions = set([r.id for r in task_model.reactions])

	# since we have multiple models, we need to evaluate the 210 tasks for each model (54 samples)
	# first, we create a structure to hold all of these results - a dictionary
	task_eval_results = {}

	# for each k (tuple with algorithm and sample information) and result (dict with reaction presences)...
	for k, result in fastcore_res_dict.items():
		# using with statements to change the COBRA model temporarily
		# this is done to knock-out reaction not appearing the FASTCORE result
		with task_model as context_specific_model:
			protected = set([k for k, v in result.items() if v])  # get reactions included in the sample-specific model
			to_remove = all_reactions - protected  # get reactions except the protected ones
			for rid in to_remove:
				context_specific_model.reactions.get_by_id(rid).knock_out()  # knock-out reactions not in the model

			# create a task evaluator instance with the context specific model and the supplied task list and solver
			task_eval = TaskEvaluator(model=context_specific_model, tasks=task_list, solver='CPLEX')

			# get task names (for future reference)
			task_names = task_eval.tasks

			# use the batch_function from the TaskEvaluator class (takes the name of a loaded task, a params
			# dictionary with the task evaluator associated to the 'tev' key) and set the amount of threads to be used
			batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names, {'tev': task_eval}, threads=10)
		# each element in the list of results in batch_res_tasks is a tuple of length 3 with the following:
		# 0 - boolean flag representing the task evaluation
		# 1 - Solution instance used to evaluate the task
		# 2 - A dictionary with reactions supposed to be active mapped to True/False according to that criterion

		# keep only items 0 and 2 of the task result - we don't need the flux distribution
		task_csm_res = {k: (v[0], v[2]) for k, v in dict(zip(task_names, batch_res_tasks)).items()}
		# assign this dictionary to it's sample on the master results dictionary
		task_eval_results[k] = task_csm_res

	# save these results for later analysis as a JSON file
	with open(TASK_RESULTS_PATH, 'w') as f:
		f.write(JSONEncoder().encode(task_eval_results))
