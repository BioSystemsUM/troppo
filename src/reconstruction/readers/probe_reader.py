"""
 Created by Jorge Gomes on 19/03/2018
 source
 probe_reader
 
"""
import numpy as np
from pandas import read_csv


class ProbeReader:
    """
    Reads expression files sourced from microarrays DBs such as Gene Expression Barcode or Gene Expression OmniBus.
    Considers each value is identified by a probeID on the first column of the file. An annotation file supplied by
    the microarray chip vendor must be supplied for appropriate probe to gene Id conversion.
    Cases where a probe has no match with convTarget nomenclature will be ignored.
    Handles cases where more than one probe translate to the same gene, and where a probe translates to more than a gene
    """

    def __init__(self, fPath, expCol, annotFile, convTarget, convSep=',', expSep=','):
        """
        Args:
            fPath: string, complete path to the file from which expresion data is read.
            expCol: int, index of the column where expression values are retrieved from.
            annotFile: string, complete path to the annotation file.
            convTarget: string, exact match to the column name of the nomenclature used for probeID to geneID conversion
                        recommended: Either Gene Symbol or Entrez Gene or equivalent.
            convSep: string, field separator used in the annotation file. Default is ",".
            expSep: string, field separator used in the probe intesity/expression file. Default is ",".
        """

        self._fpath = fPath
        self._expCol = expCol
        self._cPath = annotFile
        self._convTarget = convTarget
        self._convSep = convSep
        self._expSep = expSep
        self._IdMapping = self.__createMapping()

    def load(self):
        """
        Executes the loading of supplied omics file.

        Returns: a dictionary of geneID: expressionValue
        """
        # avoid loading when mapping does not exist
        if self._IdMapping is None:
            return

        tup_list = []  # auxiliary structure
        values = {}  # {ID: Exp_val}

        with open(self._fpath, 'r') as f:
            header = f.readline().split(self._expSep)

            if self._expCol == 0 or self._expCol > len(header):
                raise Exception('Column \'{0}\' exceeds number of columns in file, or is the probe Id column \n '
                                'Please input a valid column. File header: {1}'.format(self._expCol, header))

            else:
                for line in f:
                    fields = line.replace('\"', '').split(',')
                    genes = self._IdMapping[fields[0]].replace(' ', '').split('///')  # 1 probe : many genes
                    for geneID in genes:
                        if geneID not in ('---', ''):  # filters cases where one probe does not match an id
                            tup_list.append((geneID, float(fields[self._expCol])))
        for gene, val in tup_list:  # 2 probes translated to same gene
            if gene not in values:
                values[gene] = [val]  # simple entries
            else:
                values[gene].append(val)  # multiple entries -> mean will occur after
        new_values = {g: np.mean(val) for g, val in values.items()}
        return new_values

    # handles annotation file
    def __createMapping(self):
        field_sep = self._convSep
        mapping = {}  # handling more of one probe for the same gene

        # find header of annot file
        with open(self._cPath, 'r') as f:
            header_start = 0
            for line in f:
                if len(line.split(field_sep)) < 10:
                    header_start += 1
                else:
                    break

        # actually read the file
        annot = read_csv(self._cPath, header=header_start,sep=field_sep)
        if self._convTarget not in list(annot):
            print('convTarget is not present in the annotation file please input one of the following:','\n'
                  ,list(annot))
            return
        else:
            return dict(zip(annot.iloc[:, 0], annot[self._convTarget]))


if __name__ == '__main__':
    path = "C:/Users/Tese_Avoid_Namespaces/Tese/TsmRec/files/abc-tis-gpl570-formatted_v3.csv"
    tissue = 'brain'
    convFile = "C:/Users/Tese_Avoid_Namespaces/Tese/TsmRec/files/rembrandt_study/HG-U133_Plus_2.na35.annot.csv"
    convS = ('Probe Set ID', 'Gene Symbol', ',')
    help(ProbeReader)
    gr = ProbeReader(path, 3, convFile, convTarget="Gene Symbol")
    gr.load()


