import numpy as np
import pandas
import pandas as pd
from numpy import log


class GeneLevelThresholding:
    """
    This class is used to transform the dataframe containing the omics data and perform gene-level thresholding on omics
    data. It currently supports Global and Local thresholding approaches described by Richelle, Joshi and Lewis (2019)
    (https://doi.org/10.1371/journal.pcbi.1007185). These include:
    - global: genes with a value lower than the upper global threshold (GTU) are considered inactive; genes with a value
    greater than the lower global threshold (GTL) are considered active.
    - local t1: genes with a value lower than the upper global threshold (GTU) are considered inactive; for genes
    with a value greater than the GTU, if the value is lower than the local threshold (LT), the gene is considered
    inactive, otherwise it is considered active.
    - local t2: genes with a value lower than the upper global threshold (GTU) are considered inactive; genes with
    a value greater than the lower global threshold (GTU) are considered active; for genes with a value between the GTU
    and the lower global threshold (GTL), they are only considered active if their value is greater than the local
    threshold (LT).
    Thresholds are selected in accordance with the distribution of the data. The numbers in the thresholding options
    represent the position of the value to use. Currently, the options are: [0.1, 0.25, 0.5, 0.75, 0.9];
    the threshold value will then be the value on the dataset that corresponds to that quantile.

    Parameters
    ----------
    omics_dataframe: pandas.DataFrame
        Omics data to be thresholded.
    thresholding_strat: str
        Thresholding strategy to be used. Must be one of: global, local t1, local t2.
    global_threshold_lower: int or None, default = None
        Position of the Global Lower threshold value on the quantile list.
    global_threshold_upper: int or None, default = None
        Position of the Global Upper threshold value on the quantile list.
    local_threshold: int or None, default = None
        Position of the Local threshold value on the quantile list.
    """

    def __init__(self, omics_dataframe: pandas.DataFrame, thresholding_strat: str = 'global',
                 global_threshold_lower: int = None, global_threshold_upper: int = None, local_threshold: int = None):
        self.omics_dataframe = omics_dataframe

        self.filter_name = '_'.join(map(str, [thresholding_strat, global_threshold_lower,
                                              global_threshold_upper, local_threshold]))

        self.filter_name = self.filter_name.replace(' ', '_')

        if thresholding_strat not in ['global', 'local t1', 'local t2']:
            raise ValueError('Invalid thresholding strategy. Must be one of: global, local t1, local t2')
        self.thresholding_strat = thresholding_strat

        self.max_expression = np.log(omics_dataframe.max().max())

        qvalues = [0.1, 0.25, 0.5, 0.75, 0.9]
        self.quantiles = omics_dataframe.quantile(qvalues)
        self.global_quantiles = self.quantiles.T.apply(lambda x: x.mean())

        if thresholding_strat == 'global':
            self.global_threshold_lower = self.global_quantiles.iloc[global_threshold_lower]
            self.global_threshold_upper = None
            self.local_threshold = None
        elif thresholding_strat == 'local t1':
            self.global_threshold_lower = self.global_quantiles.iloc[global_threshold_lower]
            self.global_threshold_upper = None
            self.local_threshold = self.quantiles.iloc[local_threshold,]
        elif thresholding_strat == 'local t2':
            self.global_threshold_lower = self.global_quantiles.iloc[global_threshold_lower]
            self.global_threshold_upper = self.global_quantiles.iloc[global_threshold_upper]
            self.local_threshold = self.quantiles.iloc[local_threshold,]

    @staticmethod
    def global_thresholding(sample_series: pd.Series, gtlow: float, maxexp: float) -> dict:
        """
        Global thresholding strategy for the omics data. Processes a single sample at the time.

        Parameters
        ----------
        sample_series: pandas.Series
            Omics data from a specific sample.
        gtlow: float
            Global threshold lower value.
        maxexp: float
            Maximum expression value of the dataset.

        Returns
        -------
        filtered_sample: dict

        """
        return (sample_series / gtlow).apply(log).clip(-maxexp, maxexp).to_dict()

    @staticmethod
    def local_t1_thresholding(sample_series: pd.Series, gtlow: float, lt: pd.Series, maxexp: float) -> dict:
        """
        Local T1 thresholding strategy for the omics data. Processes a single sample at the time.

        Parameters
        ----------
        sample_series: pandas.Series
            Omics data from a specific sample.
        gtlow: float
            Global threshold lower value.
        lt: pd.Series
            Local threshold value for each sample.
        maxexp: float
            Maximum expression value of the dataset.

        Returns
        -------
        filtered_sample: dict

        """
        activity = (sample_series / gtlow).apply(log).clip(-maxexp, maxexp)
        gt_active = activity >= 0
        activity[gt_active] = (log(sample_series[gt_active] / lt[gt_active]) * maxexp).clip(-maxexp, maxexp)
        return activity.to_dict()

    @staticmethod
    def local_t2_thresholding(sample_series: pd.Series, gtlow: float, gtup: float, lt: pd.Series,
                              maxexp: float) -> dict:
        """
        Local T2 thresholding strategy for the omics data. Processes a single sample at the time.

        Parameters
        ----------
        sample_series: pandas.Series
            Omics data from a specific sample.
        gtlow: float
            Global threshold lower value.
        gtup: float
            Global threshold upper value.
        lt: pandas.Series
            Local threshold value for each gene.
        maxexp: float
            Maximum expression value of the dataset.

        Returns
        -------
        filtered_sample: dict

        """
        upp_activity = (1 + (sample_series / gtlow).apply(log)).clip(-maxexp, 1 + maxexp)
        gtu_inactive = upp_activity < 1
        low_activity = (sample_series / gtup).apply(log).clip(-maxexp, maxexp)
        gtl_maybes, gtl_lows = (low_activity >= 0) & gtu_inactive, low_activity < 0
        upp_activity[gtl_lows] = low_activity[gtl_lows]
        activity_maybe = (sample_series[gtl_maybes] / lt[gtl_maybes]). \
            apply(log).clip(-maxexp, maxexp)
        upp_activity[gtl_maybes] = activity_maybe.clip(-1, 1)
        return upp_activity.to_dict()

    def threshold_strategy(self, sample_series) -> dict:
        """
        Thresholding strategy for the omics data. Processes a single sample at the time.

        Parameters
        ----------
        sample_series: pandas.Series
            Omics data from a specific sample.

        Returns
        -------
        filtered_sample: dict

            Filtered omics data from a specific sample.

        """
        if self.thresholding_strat == 'global':
            return self.global_thresholding(sample_series, self.global_threshold_lower, self.max_expression)

        elif self.thresholding_strat == 'local t1':
            return self.local_t1_thresholding(sample_series, self.global_threshold_lower, self.local_threshold,
                                              self.max_expression)

        elif self.thresholding_strat == 'local t2':
            return self.local_t2_thresholding(sample_series, self.global_threshold_lower, self.global_threshold_lower,
                                              self.local_threshold, self.max_expression)

    def apply_thresholding_filter(self) -> pandas.DataFrame:
        """
        Thresholding filter for the omics data.

        Returns
        -------
        filtered_dataset: pandas.DataFrame
            Filtered omics dataframe.

        """
        filter_results = {}

        for sample_id in self.omics_dataframe.index:
            sample = self.omics_dataframe.loc[sample_id, :]

            filter_results[sample.name + '_' + self.filter_name] = self.threshold_strategy(sample)

        filtered_dataset = pd.DataFrame(filter_results).T

        return filtered_dataset
