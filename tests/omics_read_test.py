from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from cobra.io import read_sbml_model
from urllib.request import urlretrieve
from numpy import log, mean

if __name__ == '__main__':

	path, msg = urlretrieve('http://bigg.ucsd.edu/static/models/e_coli_core.xml')
	model = read_sbml_model(path)

	def normalize_function(df):
		dft = df.loc[:,~(df == 0).all()]
		return 5 * log(1+(dft/mean(dft, axis=0)).fillna(0))

	rdr = TabularReader('../troppo/tests/GSE135867_NQ_log_tpm_GEO.csv',ignore_samples=('gene_name'), nomenclature='any',
				omics_type='transcriptomics', dsapply=normalize_function, sample_in_rows=False, encoding='windows-1252')
	ocs = rdr.to_containers()

	rw = ReconstructionWrapper(model)
	result = rw.run_from_omics(ocs[0], 'fastcore', ('threshold', float(2*log(2))), solver='CPLEX')