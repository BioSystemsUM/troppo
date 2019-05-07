"""
 Created by Jorge Gomes on 13/03/2018
 source
 OmicsContainer
 
"""
from troppo.omics.id_converter import idConverter, searchNomenclature
import numpy as np
import copy
import re


class OmicsContainer:
    """
    OmicsContainer class to be used for the creation of objects that store omics data and other useful information, such
    as its type, and the tissue condition from where this data was obtained.
    To successfully create an OmicsContainer object one must:
        a) create an OmicsContainer object providing: a) its omictype b) the tissue/patient condition
        b) Use its .load() method providing a previously created reader object (HpaReader, ProbeReader, GenericReader)

    Once created this object can be transformed in several ways:
        a) Id conversion
        b) Value conversion
        c) Filtering by id, regular expressions, or values threshold
        d) Log transformation, or data normalization

    Main attribute is .data() which is a dictionary containing : {gene_id: Expression Value}
    """

    def __init__(self, omicstype=None, condition=None):
        """
        Creates an empty OmicsContainer object.

        Args:
            omicstype: string, defines the type of omics data present (proteomics, transcriptomics, etç)
            condition: string, contains the information describing the sample (cell type, tissue, disease)
        """
        self.otype = omicstype
        self.condition = condition
        self.data = {}
        self.nomenclature = None

    def load(self, reader):
        self.data = reader.load()
        self.nomenclature = searchNomenclature(list(self.data.keys()))

    def convertValues(self, mapping):
        """
        Converts the values in the exp_val field to different values based on a valid user supplied mapping.
        IMPORTANT: Will not work if _values contains NAs
        Mapping shall be a dictionary of either:
            - old value (may it be string or numeric): new value (may it be string or numeric)
            - tuple of (lower bound, upper bound) of old value: new value (numeric, string)
        :param mapping: a dictionary containing the mapping between the values to be converted and the desired values
        """
        if self._isNumeric():

            # range to numeric/text
            if type(list(mapping.keys())[0]) is tuple:
                if self._mapIsValid('r2n', mapping):
                    new_values = {}
                    for k, v in self.data.items():
                        for tup in mapping.keys():
                            if float(tup[0]) <= v <= float(tup[1]):
                                new_values[k] = mapping[tup]
                    self.data = new_values
                else:
                    raise Exception('Supplied mapping is not valid for the intended conversion')

            # numeric to numeric/discrete
            elif type(list(mapping.keys())[0]) in (int, float):

                # n2n
                if type(list(mapping.values())[0]) in (int, float):
                    if self._mapIsValid('n2n', mapping):
                        new_values = {k: float(mapping[v]) for k, v in self.data.items()}
                        self.data = new_values
                    else:
                        raise Exception('Supplied mapping is not valid for the intended conversion')

                # n2d
                elif type(list(mapping.values())[0]) is str:
                    if self._mapIsValid('n2d', mapping):
                        new_values = {k: mapping[v] for k, v in self.data.items()}
                        self.data = new_values
                    else:
                        raise Exception('Supplied mapping is not valid for the intended conversion')

            # avoid cases where the mapping does not follow {oldval:newval}
            else:
                raise Exception('Supplied mapping is not valid. Please supply a valid mapping with {oldval:newval}')

        # discrete to integer/float
        else:
            if self._mapIsValid('d2n', mapping):
                new_map = copy.deepcopy(mapping)  # avoid a not intended effect where this would change the user mapping

                new_map[np.NaN] = np.NaN
                new_values = {k: float(new_map[v]) for k, v in self.data.items()}
                self.data = new_values
            else:
                raise Exception('Supplied mapping is not valid for the intended conversion')

        print('Value conversion is complete!')

    def convertIds(self, new):
        """
        Redefines the ids(keys) on the data attribute.

        Args:
            new:string, designation of the new id according to hgnc
        """
        new_data = {}
        if self.nomenclature is None:
            print('No valid nomenclature was found. Please reload your data carefully')
            return
        if idConverter(self.data.keys(), self.nomenclature, new) is None:
            return
        else:
            for old, new in idConverter(self.data.keys(), self.nomenclature, new).items():
                new_data[new] = self.data[old]
            lost = len(self.data) - len(new_data.keys())

            print('ID conversion is complete! {0} entries were lost due to inexistent match in the HGNC platform'.format(lost))
            self.nomenclature = new
            self.set_data(new_data)

    def dropNA(self):  # irreversible once done
        """
        Removes every entry whose exp_val is NA
        """
        for k, v in self.data.items():
            if np.isnan(v):
                del self.data[k]

    def filterByValue(self, op, threshold):
        """
        Filters the _values attribute to match a user defined filter
        above and under use < and > operators, while between uses <= and >=.

        :param op: string, one of (above, under, between, oneof)
        :param threshold: numeric threshold for above and under, tuple of (lowerbound, upperbound) for between, string
        for included discrete levels for levels operation
        :return OmicsContainer object: a new OmicsContainer object is returned once this filter is applied. Original
        instance remains unchanged.
        """
        new_values = copy.deepcopy(self.data)

        if self._isNumeric():

            if op.lower() == 'above':
                new_values = {k: v for k, v in new_values.items() if v > threshold}

            elif op.lower() == 'under':
                new_values = {k: v for k, v in new_values.items() if v < threshold}

            elif op.lower() == 'between':
                try:
                    new_values = {k: v for k, v in new_values.items() if threshold[0] <= v <= threshold[1]}
                except IndexError:
                    print('Threshold for between operation must be a tuple of two elements')
            else:
                print('Please input a valid operation for numeric filtering: \'above\', \'under\' or \'between\'')
        else:
            if op.lower() == 'levels':
                new_values = {k: v for k, v in new_values.items() if v in threshold}
            else:
                print('Discrete filtering only supports \'levels\' operation')

        return self.__createNew(new_values)

    def filterById(self, regex):
        """
        Filters the data attribute to contain genes that match a regular expression or string supplied by the user
        :param regex: string, regular expression or string to be contained in the Gene Symbol field of the data attr.
        :return OmicsContainer object: a new OmicsContainer object is returned once this filter is applied. Original
        instance remains unchanged.
        """
        new_values = copy.deepcopy(self.data)

        try:
            exp = re.compile(regex, re.IGNORECASE)
            new_values = {k: v for k, v in new_values.items() if exp.search(k) is not None}
        except TypeError:
            print('Regex must be a string')
        return self.__createNew(new_values)

    def transform(self, func):
        """
        Applies the func to the exp_values of the data attr.
        Only compatible with numerical container.
        :param func: string, a function to be applied to the values of the container, either 'norm' or 'logx'

        Original number = x
        Transformed number x'=log10(x)
        """
        try:
            logs = {'log': np.log, 'log2': np.log2, 'log10': np.log10}

            if func.lower() in logs:
                self.data = {k: logs[func](v) for k, v in self.data.items()}

            elif func.lower() == 'norm':
                vals = [x for x in self.data.values()]
                maxV, minV = max(vals), min(vals)
                diff = maxV - minV
                self.data = {k: (v-minV)/diff for k, v in self.data.items()}
        except TypeError or ValueError:
            print('Convert to numeric values before applying normalization or log transformations')

    def _isNumeric(self):
        return type(list(self.data.values())[0]) in (int, float)

    def _mapIsValid(self, task, mapping):
        """Checks if a supllied mapping is valid, namely if all fields are present (case sensitive), and if all values
        are numerical"""
        unique = set(self.data.values())

        if task in ['n2n', 'd2n']:
            return set(mapping.keys()) == unique and len([x for x in mapping.values() if
                                                         type(x) in [int, float]]) == len(unique)
        elif task == 'n2d':
            return set(mapping.keys()) == unique and len([x for x in mapping.values() if
                                                         type(x) is str]) == len(unique)

        elif task == 'r2n':
            if len([x for x in mapping.values() if type(x) in (float, int, str)]) == len(mapping.values()):
                for v in unique:
                    covered = False
                    for tup in mapping.keys():
                        if float(tup[0]) <= v <= float(tup[1]):
                            covered = True
                    if not covered:
                        return False

                return True
            else:
                return False

    def __createNew(self, new_values):
        newOC = OmicsContainer(omicstype=self.get_OmicsType(), condition=self.get_Condition())
        newOC.set_data(new_values)
        return newOC

    # setters
    def set_type(self, newType):
        self.otype = newType

    def set_condition(self, newCond):
        self.condition = newCond

    def set_data(self, newData):
        self.data = newData

    # getters
    def __len__(self):
        return len(self.data.items())

    def get_OmicsType(self):
        return self.otype

    def get_Condition(self):
        return self.condition

    def get_Data(self):
        return self.data

    def get_Nomenclature(self):
        return self.nomenclature

    def print_values(self):
        print('Gene Symbol >>> Exp Value')
        for k, v in self.data.items():
            print("{0} >>> {1}".format(k, v))

    def __str__(self):
        return str('{0} container\n'
                   'Condition: {1}\n'
                   'Nomenclature: {2}\n'
                   '{3} Expression Values'.format(self.otype, self.condition, self.nomenclature, len(self.data.keys())))


if __name__ == '__main__':
    from troppo.readers.hpa_reader import HpaReader

    path = "../../../tests/pathology.tsv"

    d2num = {'High': 20.0,
             'Medium': 15.0,
             'Low': 10.0,
             'Not detected': -8.0}

    num2d = {20.0: 'High',
             15.0: 'Medium',
             10.0: 'Low',
             -8.0: 'Not detected'}

    n2num = {20.0: 12.0,
             15.0: 9.0,
             10.0: 5.0,
             -8.0: -1.0}

    r2n = {(0, 20): 'Positive',
           (-10, -1): 'Negative'}

    r = HpaReader(fpath=path, tissue='ovarian cancer', id_col=1, includeNA=False)
    oc = OmicsContainer(omicstype='proteomics', condition='Ovarian Cancer')
    oc.load(r)
    print(oc.get_Data())
    print(oc)

    # oc.convertValues(d2num)
    # print(oc.data)
    #
    # oc.convertValues(num2d)
    # print(oc.data)

    # ids = ["TSPAN6", "TNMD", "DPM1", "SCYL3", "C1orf112"]
    # oc.convertIds('symbol', 'entrez_id')
    # print(oc.data)
    # oc.transform('norm')
    # print(oc.values())

    # oc2 = oc.filterById('TSPAN6')
    # print(oc2.values())

    # oc2.convertIds('files/map.txt',',')
    # print(oc2.values())

    # teste values conversion
    # print(d2num)

    # oc.convertValues(d2num)
    # print(d2num)
    # print(oc.values())

    # oc.convertValues(r2n)
    # print(oc.values())

    # teste filtro valores

    #oc3 = oc.filterByValue('levels', ('Medium'))
    #print(oc3.data())

    #path = 'files/abc-tis-gpl570-formatted_v3.csv'
    #tissue = 'brain'
    #convFile = 'files/HG-U133_Plus_2.na35.annot.csv'
    #convS = ('Probe Set ID', 'Gene Symbol', ';')
    #g = GEB_Reader(path, tissue, convFile, convS)
    #oc = OmicsContainer(reader=g, Type='Transcriptomics', condition='Brain Healthy')
    #oc.print_values()
