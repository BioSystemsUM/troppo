# import numpy as np
# from numpy import min, max, isnan, mean, std
import pandas as pd


def normalize(data: pd.DataFrame) -> pd.DataFrame:
	"""
	Normalize the given data

	Parameters
	----------
	data: pd.DataFrame

	Returns
	-------
	pd.DataFrame

	"""
	# if single:
	# 	normalized_data = np.array([(n-min(data))/(max(data)-min(data)) if isnan(n) != True else -1 for n in data])
	#
	# else:
	# 	normalized_data = None

	return (data - data.min()) / (data.max() - data.min())


def z_score(data: pd.DataFrame, single=True) -> pd.DataFrame:
	"""
	Compute the Z-score of the given data

	Parameters
	----------
	data: pd.DataFrame

	single: bool, optional

	Returns
	-------
	pd.DataFrame

	"""
	# TODO this is not Z-score, is mean/std
	# if single:
	# 	z_score = np.array([(n - mean(data)) / (std(data)) if isnan(n) != True else -1 for n in data])
	# else:
	# 	z_score = None
	return (data - data.mean()) / data.std()
