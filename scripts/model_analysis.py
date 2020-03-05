from json import JSONDecoder
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy as np
from troppo.tasks.task_io import JSONTaskIO
import matplotlib.pyplot as plt

if __name__ == '__main__':
	import os
	## paths to necessary files
	## this is still hardcoded
	ROOT_FOLDER = '../troppo/scripts/ccle_breast_cancer'
	SAMPLE_INFO_PATH = os.path.join(ROOT_FOLDER,'sample_info.csv')
	MODEL_PATH = os.path.join(ROOT_FOLDER, 'Recon3D_301.xml') # SBML model path
	DATA_PATH = os.path.join(ROOT_FOLDER,'CCLE_breast_cancer_preprocessed.csv') # tabular csv with expression
	TASKS_PATH = os.path.join(ROOT_FOLDER,'nl2019_tasks_r3d_compact.json') # JSON containing metabolic tasks
	CMODL_PATH = os.path.join(ROOT_FOLDER,'Recon3D_301_consistent.xml') # Consistent SBML model path - recommended!
	TASK_RESULTS_PATH = os.path.join(ROOT_FOLDER,'results','r3d_compact_task_results_ccle_bc_new_nodrains_only_feas.json') # task evaluation
	CS_MODEL_DF_PATH = os.path.join(ROOT_FOLDER,'results','r3d_compact_ccle_bc_fastcore.csv') # context-specific models extracted from algos

	task_list = [t for t in JSONTaskIO().read_task(TASKS_PATH)]

	with open(TASK_RESULTS_PATH, 'r') as f:
		task_eval_results_slim = JSONDecoder().decode(f.read())

	fastcore_res_df = pd.read_csv(CS_MODEL_DF_PATH, index_col=[0, 1])
	fastcore_res_dict = fastcore_res_df.T.to_dict()

	k, v = zip(*task_eval_results_slim)
	task_eval_results_slim = dict(zip([tuple(i) for i in k], v))

	task_full_name_mapper = {t.name: str(t.should_fail)[0] + '-' + t.name + ' ' + t.annotations['description'] for t
							 in task_list}
	task_subsystem_mapper = {t.name: t.annotations['subsystem'] for t in task_list}

	task_res_df_orig = pd.DataFrame.from_dict(
		{k: {tk: tv[0] for tk, tv in v.items()} for k, v in task_eval_results_slim.items()}, orient='index').rename(
		columns=task_full_name_mapper)
	task_df_subsystem_series = pd.Series({fn:task_subsystem_mapper[t] for t,fn in task_full_name_mapper.items()})

	samp_info = pd.read_csv(SAMPLE_INFO_PATH, index_col=0).loc[task_res_df_orig.index.get_level_values(1),:]
	samp_info.index = pd.MultiIndex.from_arrays(zip(*[('fastcore', k) for k in samp_info.index]))

	task_res_df = task_res_df_orig.copy()
	task_res_df['subtype'] = samp_info['lineage_sub_subtype']
	df_from_dict_plot = task_res_df.groupby('subtype').sum().astype(int).apply(lambda x: x / task_res_df['subtype'].value_counts())

	df_full_plot = task_res_df_orig.T
	df_full_plot['subsystem'] = task_df_subsystem_series
	df_full_plot = df_full_plot.groupby('subsystem').mean()


	df_from_dict_subsys = df_from_dict_plot.T
	df_from_dict_subsys['subsystem'] = task_df_subsystem_series
	df_from_dict_subsys = df_from_dict_subsys.groupby('subsystem').mean()


	def cluster_heatmap(df_from_dict, plot_path, method='average', metric='jaccard', figsize=(40, 20)):
		dists_samp = dist.squareform(dist.pdist(df_from_dict.values))
		dists_task = dist.squareform(dist.pdist(df_from_dict.values.T))

		samp_linkg = sch.linkage(dists_samp, method=method, metric=metric)
		task_linkg = sch.linkage(dists_task, method=method, metric=metric)

		ord_samp = np.array(sch.dendrogram(samp_linkg)['leaves'])
		ord_task = np.array(sch.dendrogram(task_linkg)['leaves'])
		# ord_samp = sch.fcluster(samp_linkg, 0.7*max(samp_linkg[:,2]), 'distance')
		# ord_task = sch.fcluster(task_linkg, 0.7*max(task_linkg[:,2]), 'distance')
		df_from_dict_plot = df_from_dict.iloc[ord_samp, ord_task]
		plt.figure(figsize=figsize)
		plt.subplots_adjust(hspace=0.001, wspace=0.001)
		plt.gca().set_aspect(1)
		plt.pcolormesh(df_from_dict_plot)
		plt.yticks(np.arange(0.5, len(df_from_dict_plot.index), 1), [t for t in df_from_dict_plot.index])
		plt.xticks(np.arange(0.5, len(df_from_dict_plot.columns), 1), [k for k in df_from_dict_plot.columns],
				   rotation=90)
		plt.savefig(plot_path)

	cluster_heatmap(task_res_df_orig, os.path.join(ROOT_FOLDER,'results','CCLE_bc_fastcore_orig.png'),
					method='average', metric='euclidean', figsize=(30,20))

	cluster_heatmap(df_from_dict_plot, os.path.join(ROOT_FOLDER,'results','CCLE_bc_fastcore_tasks.png'),
					method='average', metric='euclidean', figsize=(30,20))

	cluster_heatmap(df_from_dict_subsys, os.path.join(ROOT_FOLDER,'results','CCLE_bc_fastcore_subsys.png'),
					method='average', metric='euclidean', figsize=(10,15))

	cluster_heatmap(df_full_plot, os.path.join(ROOT_FOLDER, 'results', 'CCLE_bc_fastcore_full_subsys.png'),
	                method='average', metric='euclidean', figsize=(15, 15))

# completed_tasks_per_model = task_res_dict.sum(axis=1)
	# fastcore_res_df.sum(axis=0).sort_values(ascending=False)
	# fastcore_res_df.sum(axis=1).sort_values()
