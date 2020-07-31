import pandas as pd
import numpy as np
from src.utilities.task import Task

class Nathan2019ConsensusPaper():

	'''
	This class has the purpose to import and make accessible the information present in the paper and in the
	supplementary material
	LINK : https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006867
	'''

	def __init__(self):
		self.root_path = 'projects/breast_mcf7/data/nl_consensus/'

	def load_medium_constraints(self):
		original = pd.read_excel(self.root_path + 'pcbi.1006867.s003.xlsx', index_col=0)
		cell_lines = [col for col in original.columns if 'Unnamed' not in col]
		gen_names = lambda x: (str(x) + '_LB', str(x) + '_UB')
		new_cols = [bound for pair in [gen_names(cell) for cell in cell_lines] for bound in pair]
		original.columns = new_cols
		final = original.drop(np.nan)
		return final

	def load_biomass(self):
		'''
		:return: biomass equations for iHsa and Recon2.2 models (based on experimental data) and cell lines growth rate
		'''
		original = pd.read_excel(self.root_path + 'pcbi.1006867.s004.xlsx', header=1)
		ihsa, recon2_2, real_biomass = original.iloc[:, [0, 1]].dropna(), original.iloc[:, [3, 4]], original.iloc[:,
																									[6, 7]].dropna()
		return ihsa, recon2_2, real_biomass

	def load_metabolic_tasks(self):
		original = pd.read_excel(self.root_path + 'pcbi.1006867.s005.xlsx')
		idx = np.where(original.ID.isna()==False)[0]
		idx = np.append(idx, original.shape[0]-1)
		return {i+1: Task(original.iloc[idx[i]:idx[i + 1], :]) for i in range(len(idx) - 1)}

	def load_hallmarks_genes(self):
		original = pd.read_excel(self.root_path + 'pcbi.1006867.s013.xlsx', header=1)
		return original

if __name__ == '__main__':
	path = './data/Nathan2019ConsensusPaper/'

	n = Nathan2019ConsensusPaper()
	x = n.load_medium_constraints()
	c, v, b = n.load_biomass()
	# z = pd.read_excel(path + 'pcbi.1006867.s005.xlsx')
	# original = pd.read_excel(path + 'pcbi.1006867.s005.xlsx')
	# idx = np.where(original.ID.isna() == False)[0]
	# idx = np.append(idx, original.shape[0] - 1)
	# tasks = {i+1: Task(original.iloc[idx[i]:idx[i + 1], :]) for i in range(len(idx) - 1)}
	original = pd.read_excel(path + 'pcbi.1006867.s013.xlsx', header=1)