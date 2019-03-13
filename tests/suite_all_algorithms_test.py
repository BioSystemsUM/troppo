import cobra
import framed
import cobamp
import pandas as pd
import numpy as np
import pickle
import scipy as sci

# for testing the algorithms
from cobamp.wrappers import COBRAModelObjectReader
from reconstruction.methods.tINIT import tINIT
from reconstruction.methods.fastcore import FASTcore
from reconstruction.reconstruction_properties import FastcoreProperties, tINITProperties

if __name__ == '__main__':
	from random import random

	rootpath = 'E:/framed_hm/sbml_files/'
	testpath = 'D:/Matlab/'

	def sample_models():
		# # model1 = cobra.io.read_sbml_model(rootpath + 'Recon2_mod_gpr.xml')
		# # model1 = cobra.io.read_sbml_model(rootpath + 'Recon2_v04.xml')
		model1 = cobra.io.read_sbml_model('tests/test_model_fastcore.xml')
		# # # model2 = framed.load_cbmodel(rootpath+'Recon2_mod_gpr.xml')
		# #
		# # with open('recon2.pkl','wb') as f:
		# # 	pickle.dump(model1, f)
		#
		# # with open('E:/reconstruction/recon2.pkl', 'rb') as f:
		# # 	model1 = pickle.load(f)
		#
		reader1 = cobamp.wrappers.COBRAModelObjectReader(model1)
		#
		reader1.get_irreversibilities(True)
		#
		S = reader1.get_stoichiometric_matrix()
		lb, ub = reader1.get_model_bounds(False, True)
		rx_names = reader1.get_metabolite_and_reactions_ids()[0]
		#
		#
		# # reader2 = cobamp.wrappers.FramedModelObjectReader(model2)
		#
		# # model2.reactions['R_ATPS4m'].gpr.to_string()
		#
		# def test_reader(reader):
		# 	reader.get_model_gprs()
		# 	reader.convert_gprs_to_list('ATPS4m')
		# 	expression = {k: random() * 10 for k in reader.get_model_genes()}
		# 	reader.g2rx(expression, min, max)
		#
		# Generate random Data
		# reader1.get_model_gprs()[123]
		expression = {k: random() * 10 for k in reader1.get_model_genes()}
		# reader1.convert_gprs_to_list(reader1.get_rx_instances()[123].id)
		map = reader1.g2rx(expression, min, max)
		map_df = pd.DataFrame.from_dict(map, orient='index')
		# pd.Series(map_df[map_df[0] > 5].index.values).to_csv(testpath + 'reac_to_test.csv', index=False, header=True)
		# map_df = pd.Series(map_df[map_df[0] > 5].index.values)
		core_idx = [reader1.reaction_id_to_index(x) for x in map_df[map_df[0] > 5].index.values]

		return model1, reader1, S, lb, ub, rx_names, map, map_df, core_idx


	def consistent_recon2_from_mat():
		x = sci.io.loadmat('D:/Matlab/cobratoolbox/test/verifiedTests/analysis/testFASTCORE/FastCoreTest.mat')
		model = x['ConsistentRecon2']
		S = model['S'][0][0].toarray()
		lb = model['lb'][0][0].ravel()
		ub = model['ub'][0][0].ravel()
		# reactions = model['rxns'][0][0].ravel()
		reactions = [x[0][0] for x in model['rxns'][0][0]]
		reac_scores = np.array([random() * 10 for k in reactions])
		# map = {k: random() * 10 for k in reactions}
		# map_df = pd.DataFrame.from_dict(map, orient='index')
		# map_idx = {k:map[k] for k in reactions}

		return S, lb, ub, reactions, reac_scores

	def testFastcore():
		matlab_core = pd.read_csv('./tests/fastcore_core_test_matlab.csv', header=None) - 1

		x = sci.io.loadmat('D:/Matlab/cobratoolbox/test/verifiedTests/analysis/testFASTCORE/FastCoreTest.mat')
		model = x['ConsistentRecon2']
		S = model['S'][0][0].toarray()
		lb = model['lb'][0][0].ravel()
		ub = model['ub'][0][0].ravel()

		## Testing the algorithms with the matlab versions
		# FASTCORE
		print('testing FASTCORE')
		f = FASTcore(S, lb, ub, FastcoreProperties(solver = 'GUROBI', core=list(matlab_core[0].values)))

		tissue_reactions = f.fastcore()
		print(len(tissue_reactions), tissue_reactions)
		# irrev_matlab = pd.read_csv('E:/reconstruction/tests/irrev_matlab.csv', header = None)
		# np.setdiff1d(tissue_reactions, irrev_matlab[0]-1)
		# pd.DataFrame.from_dict(map,orient='index').to_csv(testpath + 'map_to_test.csv')#,header=False, na_rep=0)

	# tINIT
	def test_tINIT():
		S, lb, ub, reactions, reac_scores = consistent_recon2_from_mat()
		reac_scores.T.tofile('D:/Matlab/tINIT_map.csv', sep= ',')

		print('testing tINIT')
		t = tINIT(S, np.array(lb), np.array(ub),
				  tINITProperties(reactions_scores=reac_scores, present_metabolites=[], essential_reactions=[],
								  production_weight=0.5, allow_excretion=False, no_reverse_loops=False, solver = "GUROBI"))
		t.preprocessing()
		t.build_problem()
		res = t.solve_problem()
		res.sort()
		print(np.int_(res) + 1)
		print(len(np.int_(res) + 1))


	### Tests to Run ###

	test_tINIT()
