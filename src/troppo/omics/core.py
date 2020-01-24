import copy
import re

import numpy as np

from troppo.omics.id_converter import searchNomenclature, idConverter


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

    def __init__(self, omicstype=None, condition=None, data=None, nomenclature=None):
        """
        Creates an empty OmicsContainer object.

        Args:
            omicstype: string, defines the type of omics data present (proteomics, transcriptomics, etç)
            condition: string, contains the information describing the sample (cell type, tissue, disease)
        """
        self.otype = omicstype
        self.condition = condition
        self.nomenclature = nomenclature
        if data == None:
            self.data = {}
        else:
            self.load(data)

    def load(self, arg, **kwargs):
        if isinstance(arg, dict):
            self.data = arg
        else:
            self.data = arg.load(**kwargs)
        if self.nomenclature is None:
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


    def get_integrated_data_map(self, model_reader, and_func=min, or_func=max):
        """
        Function responsible for the integration of different omics data with a metabolic model loaded with framed package.
        Matches model ids for gene_ids, metabolites or reaction ids with those present in the omicsContainer object.

        :param model_reader: (obj) a cobamp AbstractModelObjectReader object
        :param and_func:(func) the mathematical function to replace the "AND" operator present in the Gene-Protein-Rules
        :param or_func:(func) the mathematical function to replace the "OR" operator present in the Gene-Protein-Rules


        :return m: (obj) an OmicsDataMap object which contains the mapping between reactions/metabolites and its fluxes
        based on the supplied omics data.
        """

        def g2rIntegrate():
            """
            Handles integration of both proteomics and transcriptomics expression data relying on framed's gene2reaction
            """
            # suffixAndPrefix()
            d = model_reader.get_reaction_scores(self.get_Data(), or_fx=or_func, and_fx=and_func)
            return aux_createMap(d, 'ReactionDataMap')

        def aux_createMap(mMap, mapType):
            m = OmicsDataMap(mMap, mapType)
            return m

        # execution commands

        omicsType = self.otype.lower()

        if omicsType.lower() in ['proteomics', 'transcriptomics']:
            return g2rIntegrate()
        else:
            raise Exception('Omics data type not yet supported')


    def print_values(self):
        print('Gene Symbol >>> Exp Value')
        for k, v in self.data.items():
            print("{0} >>> {1}".format(k, v))

    def __str__(self):
        return str('{0} container\n'
                   'Condition: {1}\n'
                   'Nomenclature: {2}\n'
                   '{3} Expression Values'.format(self.otype, self.condition, self.nomenclature, len(self.data.keys())))


class OmicsDataMap:
    """
    Stores integrated omics data, matching a given metabolic model

    """

    def __init__(self, scores, mapType):
        self._mapType = mapType
        self._scores = scores

    # getters
    def __len__(self):
        return len(self._scores.items())

    def __str__(self):
        return self.get_scores()

    def mapType(self):
        return self._mapType

    def get_scores(self):
        return self._scores

    def select(self, op, threshold):
        """
        Filtering the original reaction scores to be under or above a threshold. Above or under operations use the
        >= and <= operators

        Args:
            op: str, either "above" or "under" determining which scores shall be chosen
            threshold: num, either a float or an integer whether under or above all scores shall be chosen


        """

        if type(threshold) not in [int, float]:
            print('Select threshold must be numeric!')
            return

        if op.lower() == 'above':
            res = {x: y for x, y in self._scores.items() if y is not None and y >= threshold}

        elif op.lower() == 'under':
            res = {x: y for x, y in self._scores.items() if y is not None and y <= threshold}

        else:
            print('Select operation must be either \'above\' or \'under\'')
            return

        # self.set_scores(res)
        return set(res.keys())

    # setters
    def set_scores(self, newScores):
        self._scores = newScores


