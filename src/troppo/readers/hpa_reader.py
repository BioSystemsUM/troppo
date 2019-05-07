"""
 Created by Jorge Gomes on 09/03/2018
 source
 HPA_Reader
 
"""
import numpy as np


class HpaReader:
    """
    Reads the HPA pathology.tsv file from a fpath in the system.
    Discrete values are converted to numerical and expression values account for the level with the most patients
    """

    def __init__(self, fpath, tissue, id_col=0, includeNA=False):
        """
        Args:
            fpath: string, complete path to the file from which omics data is read
            tissue: string, Exactly as in the file, regarding the column where expression values should be retrieved
            id_col: int, either 0 (="ensembl") or 1(="gene_symbol") regarding which column shall be used for gene id
            includeNA: boolean, flag if NA values should be included or not
        """
        self._tissue = tissue
        self._id_col = id_col
        self._path = fpath
        self._includeNA = includeNA

    def load(self):
        """
            Executes the loading of supplied omics file.
            Returns: a dictionary of geneID: expressionValue
            """
        if self._id_col not in (0, 1):
            print('Invalid id_col. Using column 0 for gene ids')
            self._id_col = 0
        with open(self._path, 'r') as f:
            header = f.readline().split('\t')
            levels = header[3:7]

            record = {}  # {Gene symbol: Expression Value}

            for line in f:
                fields = line.split('\t')
                if fields[2] == self._tissue:
                    # record
                    if not np.isnan(_handle_exp_val(fields[3:7])):
                        record[fields[self._id_col]] = levels[_handle_exp_val(fields[3:7])]
                    elif self._includeNA:
                        record[fields[self._id_col]] = np.NaN

        return record


# Auxiliary functions


def _handle_exp_val(exp_values):
    """Retrieves the index of the expression value with the most patients."""

    if exp_values == ['', '', '', '']:
        return np.NaN
    else:
        max_idx = [i for i, x in enumerate(exp_values) if x == max(exp_values)]
        return max_idx[0]


# as of now not being used
def _handle_prog(prog):
    """Retrieves the output prognostic based on the score placement in HPA file"""
    # record['Prognostic'].append(progs[handle_prog(fields[7:])].strip('\n') if handle_prog(fields[7:])
    # is not 'None' else 'None')
    if prog == ['', '', '', '\n']:
        return 'None'
    else:
        idx = [len(i) for i in prog]
        return idx.index(max(idx))


if __name__ == '__main__':
    PATH = "C:/Users/Tese_Avoid_Namespaces/Tese/TsmRec/files/pathology.tsv"
    d2num = {'High': 20.0,
             'Medium': 15.0,
             'Low': 10.0,
             'Not detected': -8.0}

    hpa = HpaReader(PATH, 'breast cancer', id_col=2, includeNA=False)
    a = hpa.load()
    print(a)


