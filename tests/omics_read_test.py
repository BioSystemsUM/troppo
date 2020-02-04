from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.tasks.core import TaskEvaluator
from troppo.tasks.task_io import JSONTaskIO
from cobra.io import read_sbml_model
from cobra.flux_analysis.variability import find_blocked_reactions
from cobamp.utilities.parallel import batch_run
import pandas as pd

from numpy import log
import re

if __name__ == '__main__':

	## helper functions
	patt = re.compile('__COBAMPGPRDOT__[0-9]{1}') # find .{number} references
	replace_alt_transcripts = lambda x: patt.sub('', x) # replace .{number} with nothing


	MODEL_PATH = '../troppo/tests/Recon3D_301.xml'
	DATA_PATH = '../troppo/tests/ccle_breast_cancer/CCLE_breast_cancer_preprocessed.csv'
	TASKS_PATH = '../troppo/tests/nl2019_tasks_r3d_compact.json'

	## model preprocessing
	model = read_sbml_model(MODEL_PATH)
	with model as m:
		for k in m.boundary:
			k.bounds = -1000, 1000
		blocked_reactions = find_blocked_reactions(m)
	model.remove_reactions(blocked_reactions, remove_orphans=True)

	for k in model.boundary:
		k.bounds = -1000, 1000


	ocs = TabularReader(path=DATA_PATH, nomenclature='entrez_id', omics_type='transcriptomics').to_containers()
	rw = ReconstructionWrapper(model, ttg_ratio=9999, gpr_gene_parse_function=replace_alt_transcripts)

	t = (5 * log(2))

	omics_sample = ocs[0]
	fastcore_result = rw.run_from_omics(omics_container=ocs[0], algorithm='fastcore',
										integration_strategy=('threshold', t), solver='CPLEX')

	x = pd.DataFrame.from_dict({('fastcore','a'): fastcore_result, ('fastcore','b'): fastcore_result}, orient='index')

	task_list = JSONTaskIO().read_task(TASKS_PATH)

	drains = set([r.id for r in model.boundary])
	protected = drains | set([k for k,v in fastcore_result.items() if v])
	to_remove = set([r.id for r in model.reactions]) - protected

	def context_apply(model):
		for k in to_remove:
			model.set_reaction_bounds(index=k, lb=0, ub=0)

	with model as context_specific_model:
		for rid in to_remove:
			context_specific_model.reactions.get_by_id(rid).knock_out()
		task_eval = TaskEvaluator(model=context_specific_model, tasks=task_list, solver='CPLEX')
		batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_eval.tasks, {'tev': task_eval}, threads=10)

	task_eval_dict = dict(zip(task_eval.tasks, batch_res_tasks))

	# slow af
	# sign_change = lambda t: lambda x: {k: float(0 if v is None else v/t if (v >= t) else -(1 - (v/t))) for k,v in x.items()}
	# init_result = rw.run_from_omics(omics_container=ocs[0], algorithm='tinit',
	# 									integration_strategy=('continuous', sign_change(t)), solver='CPLEX')


	# score_in_interval = lambda x,y: lambda dm: {k for k,z in dm.get_scores().items() if z is not None and x < z <= y}
	# corda_thresholds = [score_in_interval(*k) for k in [(t, 9999),(t/2, t),(0, t/4)]]
	# corda_result = rw.run_from_omics(omics_container=ocs[0], algorithm='corda',
	# 								 	integration_strategy=('custom', corda_thresholds))


	## CORSO test
