"""
 Created by Jorge Gomes on 06/09/2018
 source
 generic_reader
 
"""
from pandas import read_csv, DataFrame

from ..core import OmicsContainer


class TabularReader(object):
    def __init__(self, path_or_df, index_col=0, sample_in_rows=True, header_offset=0, cache_df=False, ignore_samples=None,
                 omics_type='transcriptomics', nomenclature=None, dsapply=None, **kwargs):
        self.path, self.index_col, self.sample_axis, self.header_offset = \
            path_or_df, index_col, sample_in_rows, header_offset
        self.pandas_args = kwargs
        self.dsapply = dsapply
        self.dfcache = None
        self.cache_df = cache_df
        self.ignore_samples = ignore_samples
        self.omics_type = omics_type
        self.nomenclature = nomenclature

    def __iter__(self):
        if self.dfcache is None:
            if isinstance(self.path, DataFrame):
                df = self.path
            else:
                df = read_csv(self.path, index_col=self.index_col, **self.pandas_args)
            self.dfcache = df
        elif self.cache_df:
            df = self.dfcache

        else:
            df = self.dfcache

        if not self.sample_axis:
            df = df.T

        if self.ignore_samples is not None and len(self.ignore_samples) > 0:
            df = df.drop(labels=self.ignore_samples, axis=0)

        if self.dsapply is not None:
            df = self.dsapply(df)


        for name, series in df.iterrows():
            yield (name, series.to_dict())

    def to_containers(self):
        ocs = [OmicsContainer(data=data, condition=sample, nomenclature=self.nomenclature, omicstype=self.omics_type)
               for sample, data in self]
        self.dsapply = None
        return ocs

class GenericReader:
    """
    A generic reader to be used with omics files that are unable to be loaded by ProbeReader, or HpaReader, such as
    RNA-seq files from the gdc. Capable of handling files with additional info before the file header when supplied
    header_start by the user.
    """
    def __init__(self, path, idCol, expCol, header_start=0, sep=','):
        """
        Args:
            path: string, complete path to the file from which expresion data is read.
            idCol: integer or string, either the name of the identifier column or its index in the file header
            expCol: integer or string, either the name of the expression values column or its index in the file header
            header_start: integer, line of the file header. Default = 0
            sep: string, field separator used in the omics file. Default = ","
        """
        self._path = path
        self._idCol = idCol
        self._expCol = expCol
        self._headerStart = header_start
        self._sep = sep

    def load(self, **kwargs):
        """
        Executes the loading of supplied omics file.

        Returns: a dictionary of geneID: expressionValue
        """
        omics_ds = read_csv(self._path, header=self._headerStart, sep=self._sep)
        omics = {}
        t_id = type(self._idCol)
        t_exp = type(self._expCol)

        if t_id not in (int, str) or t_exp not in (int, str):
            print('Invalid idCol or expCol. Please input one of the following:\n', omics_ds.columns.values)
            return
        try:
            if t_id == t_exp == int:
                omics = dict(zip(omics_ds.iloc[:, self._idCol], omics_ds.iloc[:, self._expCol]))
            elif t_id == t_exp == str:
                omics = dict(zip(omics_ds[self._idCol], omics_ds[self._expCol]))
            elif t_id is int and t_exp is str:
                omics = dict(zip(omics_ds.iloc[:, self._idCol], omics_ds[self._expCol]))
            elif t_id is str and t_exp is int:
                omics = dict(zip(omics_ds[self._idCol], omics_ds.iloc[:, self._expCol]))
        except (KeyError, IndexError):
            print('One or both of the supplied columns do not match any column of the file or supplied column index is'
                  ' out of range\n Header Size:{0}\n Columns:{1}'
                  .format(len(omics_ds.columns.values)-1, omics_ds.columns.values))

        return omics



if __name__ == "__main__":
    # path1 = "C:/Users/Tese_Avoid_Namespaces/Tese/TsmRec/files/abc-gpl571-formatted_v3.csv"
    # gr = GenericReader(path1,'probe_id', 22)
    # gr.load()
    reader = TabularReader('tests/small_transcriptomics_dataset.csv')

