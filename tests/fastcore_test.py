from cobra.io import read_sbml_model

from cobamp.wrappers import COBRAModelObjectReader
from troppo.methods.fastcore import FASTcore, FastcoreProperties
import pandas as pd
from urllib.request import urlretrieve

import numpy as np

if __name__ == '__main__':

	path, content = urlretrieve('http://bigg.ucsd.edu/static/models/iAF1260.xml')

	model = read_sbml_model(path)
	cobamp_model = COBRAModelObjectReader(model)

	cobamp_model.get_irreversibilities(True)

	S = cobamp_model.get_stoichiometric_matrix()
	lb, ub = cobamp_model.get_model_bounds(False, True)
	rx_names = cobamp_model.get_metabolite_and_reactions_ids()[0]

	# S = np.array([[1, -1, 0, 0, -1, 0, -1, 0, 0],
	# 			  [0, 1, -1, 0, 0, 0, 0, 0, 0],
	# 			  [0, 1, 0, 1, -1, 0, 0, 0, 0],
	# 			  [0, 0, 0, 0, 0, 1, -1, 0, 0],
	# 			  [0, 0, 0, 0, 0, 0, 1, -1, 0],
	# 			  [0, 0, 0, 0, 1, 0, 0, 1, -1]])
	# M, N = S.shape
	# lb = np.array([0] * N).astype(float)
	# lb[3] = -10
	# ub = np.array([10] * N).astype(float)
	#
	# names = ['R' + str(i + 1) for i in range(N)]

	f = FASTcore(S, lb, ub, FastcoreProperties(core=[1004]))
	# f = FASTcore(S, lb, ub, FastcoreProperties(core=[8]))
	tissue_reactions = f.fastcore()

	matlab_test = pd.read_csv('./tests/fastcore_matlab_iAF1260.txt', header=None)
	matlab_test_array = np.array(matlab_test.loc[:,0])
	matlab_test_array_idx = [cobamp_model.reaction_id_to_index(i) for i in matlab_test_array]

	diff = np.setdiff1d(tissue_reactions, matlab_test_array_idx)
	diff2 = np.setdiff1d(matlab_test_array_idx, tissue_reactions)
	[rx_names[i] for i in diff]
	[rx_names[i] for i in diff2]

	print([i + 1 for i in tissue_reactions])
	names_to_test = [rx_names[i] for i in f.properties['core_idx']]
	with model:
		model.remove_reactions([rx_names[int(i)] for i in range(S.shape[1]) if i not in tissue_reactions])
		sol1 = model.optimize()
		print(sol1.to_frame())
		print(sol1.to_frame().loc[pd.Series(names_to_test),:])

	with model:
		model.remove_reactions([rx_names[int(i)] for i in range(S.shape[1]) if i not in matlab_test_array_idx])
		sol2 = model.optimize()
		print(sol2.to_frame())
		print(sol2.to_frame().loc[pd.Series(names_to_test),:])

	sol_f = sol1.to_frame().join(sol2.to_frame(), lsuffix = str(1), rsuffix = str(2), how = 'outer')
	sol_f['diff'] = (sol_f['fluxes1'] - sol_f['fluxes2']).abs()
	sol_f[sol_f['diff'] >= 1e3]