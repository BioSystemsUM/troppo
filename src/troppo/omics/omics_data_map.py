"""
 Created by Jorge Gomes on 28/03/2018
 TsmRec
 OmicsDataMap
 
"""


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
            res = {x: y for x, y in self._scores.items() if y >= threshold}

        elif op.lower() == 'under':
            res = {x: y for x, y in self._scores.items() if y <= threshold}

        else:
            print('Select operation must be either \'above\' or \'under\'')
            return

        self.set_scores(res)
        return set(self._scores.keys())

    # setters
    def set_scores(self, newScores):
        self._scores = newScores


if __name__ == '__main__':
    import pickle
    from readers.hpa_reader import HpaReader
    from omics.omics_container import OmicsContainer
    from omics.integrate import integrateOmics
    import time

    now = time.time()
    # metabolic model
    pickle_file_simplified = "../../models/recon1_pickle_simplified"
    pickle_model = open(pickle_file_simplified, 'rb')
    model = pickle.load(pickle_model)
    print('Model loaded')
    print(time.time() - now)

    # reader
    omics_file = "C:/Users/Tese_Avoid_Namespaces/Tese/TsmRec/files/pathology.tsv"
    d2num = {'High': 20.0,
             'Medium': 15.0,
             'Low': 10.0,
             'Not detected': -8.0}
    hpa = HpaReader(omics_file, 'breast cancer', id_col=1, includeNA=False)
    print('Omics read')
    print(time.time() - now)

    # omicsContainer
    oc = OmicsContainer(omicstype='transcriptomics', condition='breast_cancer')
    oc.load(hpa)
    oc.convertValues(d2num)
    print(oc)

    print(oc.nomenclature)
    #print(len(oc.data))
    print('container created')
    print(time.time() - now)

    # integration
    dm = integrateOmics(model, oc)
    print(dm.get_scores())
    print(len(dm.get_scores()))
    print('integration done')
    print(time.time() - now)