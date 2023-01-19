import abc

MINSUM = (min, sum)
MINMAX = (min, max)

class ScoreIntegrationStrategy():
    __metaclass__ = abc.ABCMeta
    '''

    troppo.omics.core.OmicsDataMap goes in
    some array etc... comes out

    '''

    @staticmethod
    @abc.abstractmethod
    def integrate(self, data_map): pass


class ReactionProtectionMixin():
    def __init__(self, protected_reactions):
        self.protected_reactions = protected_reactions


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


class AdjustedScoreIntegrationStrategy(ScoreIntegrationStrategy, ReactionProtectionMixin):
    def __init__(self, protected_reactions):
        super().__init__(protected_reactions)

    def integrate(self, data_map):
        maxv = max([k for k in data_map.get_scores().values() if k is not None])
        scores = {k: (v / maxv if v < 0 else v) if v is not None else 0 for k, v in data_map.get_scores().items()}
        scores.update({x: max(scores.values()) for x in self.protected_reactions})
        return scores

class DefaultCoreIntegrationStrategy(ScoreIntegrationStrategy, ReactionProtectionMixin):
    def __init__(self, threshold, protected_reactions):
        super().__init__(protected_reactions)
        self.__threshold = threshold

    def integrate(self, data_map):
        return [[k for k, v in data_map.get_scores().items() if
                 (v is not None and v > self.__threshold) or k in self.protected_reactions]]

class ThresholdSelectionIntegrationStrategy(ScoreIntegrationStrategy):
    def __init__(self, thresholds):
        if isinstance(thresholds, (int, float)):
            self.thresholds = [thresholds]
        else:
            self.thresholds = thresholds

    def integrate(self, data_map):
        tvals = [data_map.select(op='above', threshold=float(t)) for t in self.thresholds]
        return tvals[0] if len(tvals) < 2 else tvals
