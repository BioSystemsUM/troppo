import abc

class ScoreIntegrationStrategy():
	__metaclass__ = abc.ABCMeta
	'''

	troppo.omics.core.OmicsDataMap goes in
	some array etc... comes out

	'''

	@staticmethod
	@abc.abstractmethod
	def integrate(self, data_map): pass


class ContinuousScoreIntegrationStrategy(ScoreIntegrationStrategy):
	def __init__(self, score_apply=None):
		self.score_apply = score_apply

	def integrate(self, data_map):
		return data_map.get_scores() if self.score_apply is None else self.score_apply(data_map.get_scores())


class CustomSelectionIntegrationStrategy(ScoreIntegrationStrategy):
	## TODO: group_functions must be a dict
	def __init__(self, group_functions):
		self.group_functions = group_functions

	def integrate(self, data_map):
		## TODO: return type must be a dict(str -> array)
		tvals = [f(data_map) for f in self.group_functions]
		return tvals[0] if len(tvals) < 2 else tvals

class ThresholdSelectionIntegrationStrategy(ScoreIntegrationStrategy):
	def __init__(self, thresholds):
		if isinstance(thresholds, (int, float)):
			self.thresholds = [thresholds]
		else:
			self.thresholds = thresholds

	def integrate(self, data_map):
		tvals = [data_map.select(op='above', threshold=float(t)) for t in self.thresholds]
		return tvals[0] if len(tvals) < 2 else tvals
