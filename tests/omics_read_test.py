from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from cobra.io import read_sbml_model
from numpy import log
import re

if __name__ == '__main__':

	## helper functions
	patt = re.compile('_AT[0-9]{1}') # find .{number} references
	replace_alt_transcripts = lambda x: patt.sub('', x) # replace .{number} with nothing


	MODEL_PATH = '../ccle_model_reconstruction/resources/models/Recon3D.xml'
	DATA_PATH = '../troppo/tests/ccle_breast_cancer/CCLE_breast_cancer_preprocessed.csv'

	model = read_sbml_model(MODEL_PATH)
	ocs = TabularReader(path=DATA_PATH, nomenclature='entrez_id', omics_type='transcriptomics').to_containers()

	rw = ReconstructionWrapper(model, ttg_ratio=9999, gpr_gene_parse_function=replace_alt_transcripts)
	t = (5 * log(2))

	#
	# fastcore_result = rw.run_from_omics(omics_container=ocs[0], algorithm='fastcore',
	# 									integration_strategy=('threshold', t), solver='CPLEX')

	# slow af
	# sign_change = lambda t: lambda x: {k: float(0 if v is None else v/t if (v >= t) else -(1 - (v/t))) for k,v in x.items()}
	# init_result = rw.run_from_omics(omics_container=ocs[0], algorithm='tinit',
	# 									integration_strategy=('continuous', sign_change(t)), solver='CPLEX')


	score_in_interval = lambda x,y: lambda dm: {k for k,z in dm.get_scores().items() if z is not None and x < z <= y}
	corda_thresholds = [score_in_interval(*k) for k in [(t, 9999),(t/2, t),(0, t/4)]]
	corda_result = rw.run_from_omics(omics_container=ocs[0], algorithm='corda',
									 	integration_strategy=('custom', corda_thresholds))


	## CORSO test
