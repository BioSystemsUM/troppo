import abc

from troppo.omics.core import OmicsDataMap

MINSUM = (min, sum)
MINMAX = (min, max)


class ScoreIntegrationStrategy:
    __metaclass__ = abc.ABCMeta
    '''
    This class is used to integrate the scores of the different omics data.
    
    Attributes
    ----------
    data_map : OmicsDataMap
        The data map containing the gene scores to be integrated into reaction scores.

    '''

    @staticmethod
    @abc.abstractmethod
    def integrate(self, data_map: OmicsDataMap): pass


class ReactionProtectionMixin:
    """
    This class is used to protect reactions from being removed by the integration strategy.

    Attributes
    ----------
    protected_reactions : list
        The list of reactions to be protected from being removed by the integration strategy.

    """
    def __init__(self, protected_reactions: list):
        self.protected_reactions = protected_reactions


class ContinuousScoreIntegrationStrategy(ScoreIntegrationStrategy):
    """
    This class is used to integrate continuous scores.

    Attributes
    ----------
    score_apply : function
        The function to be applied to the scores.

    """
    def __init__(self, score_apply=None):
        self.score_apply = score_apply

    def integrate(self, data_map: OmicsDataMap) -> dict:
        """
        This method is used to integrate the scores of the different omics data.

        Parameters
        ----------
        data_map: OmicsDataMap
            The data map containing the gene scores to be integrated into reaction scores.

        Returns
        -------
        dict: The integrated scores.
        """
        return data_map.get_scores() if self.score_apply is None else self.score_apply(data_map.get_scores())


class CustomSelectionIntegrationStrategy(ScoreIntegrationStrategy):
    """
    This class is used to integrate the scores of the different omics data.

    Attributes
    ----------
    group_functions : dict
        The dictionary containing the functions to be applied to the scores.

    """
    ## TODO: group_functions must be a dict
    def __init__(self, group_functions: dict):
        self.group_functions = group_functions

    def integrate(self, data_map: OmicsDataMap) -> dict:
        """
        This method is used to integrate the scores of the different omics data.

        Parameters
        ----------
        data_map: OmicsDataMap
            The data map containing the gene scores to be integrated into reaction scores.

        Returns
        -------
        list: The integrated scores.
        """
        ## TODO: return type must be a dict(str -> array)
        tvals = [f(data_map) for f in self.group_functions]
        return tvals[0] if len(tvals) < 2 else tvals


class AdjustedScoreIntegrationStrategy(ScoreIntegrationStrategy, ReactionProtectionMixin):
    """
    This class is used to integrate the scores of the different omics data.

    Attributes
    ----------
    protected_reactions : list
        The list of reactions to be protected from being removed by the integration strategy.

    """
    def __init__(self, protected_reactions: list):
        super().__init__(protected_reactions)

    def integrate(self, data_map: OmicsDataMap) -> dict:
        """
        This method is used to integrate the scores of the different omics data.

        Parameters
        ----------
        data_map: OmicsDataMap
            The data map containing the gene scores to be integrated into reaction scores.

        Returns
        -------
        dict: The integrated scores.
        """
        maxv = max([k for k in data_map.get_scores().values() if k is not None])
        scores = {k: (v / maxv if v < 0 else v) if v is not None else 0 for k, v in data_map.get_scores().items()}
        scores.update({x: max(scores.values()) for x in self.protected_reactions})
        return scores


class DefaultCoreIntegrationStrategy(ScoreIntegrationStrategy, ReactionProtectionMixin):
    """
    This class is used to integrate the scores of the different omics data.

    Attributes
    ----------
    threshold: float or int
        The threshold to be applied to the scores.

    protected_reactions : list
        The list of reactions to be protected from being removed by the integration strategy.

    """
    def __init__(self, threshold: float or int, protected_reactions: list):
        super().__init__(protected_reactions)
        self.__threshold = threshold

    def integrate(self, data_map: OmicsDataMap) -> list:
        return [[k for k, v in data_map.get_scores().items() if
                 (v is not None and v > self.__threshold) or k in self.protected_reactions]]


class ThresholdSelectionIntegrationStrategy(ScoreIntegrationStrategy):
    """
    This class is used to integrate the scores of the different omics data.

    Attributes
    ----------
    thresholds : list or float or int
        The thresholds to be applied to the scores. If a list is provided, the integration will be performed for each
        threshold. If a single value is provided, the integration will be performed only once.

    """
    def __init__(self, thresholds: list or float or int):
        if isinstance(thresholds, (int, float)):
            self.thresholds = [thresholds]
        else:
            self.thresholds = thresholds

    def integrate(self, data_map: OmicsDataMap) -> list:
        tvals = [data_map.select(op='above', threshold=float(t)) for t in self.thresholds]
        return tvals[0] if len(tvals) < 2 else tvals
