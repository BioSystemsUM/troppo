from cobra.io import read_sbml_model

from cobamp.wrappers.cobra import COBRAModelObjectReader
from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties
import pandas as pd
from urllib.request import urlretrieve

if __name__ == '__main__':
    path, content = urlretrieve('http://bigg.ucsd.edu/static/models/iAF1260.xml')

    model = read_sbml_model(path)
    cobamp_model = COBRAModelObjectReader(model)

    cobamp_model.get_irreversibilities(True)

    s_matrix = cobamp_model.get_stoichiometric_matrix()
    lower_bound, upper_bound = cobamp_model.get_model_bounds(False, True)
    reaction_names = cobamp_model.get_reaction_and_metabolite_ids()[0]

    fastcore_algorithm = FASTcore(s_matrix, lower_bound, upper_bound, FastcoreProperties(core=[8], solver='CPLEX'))
    tissue_reactions = fastcore_algorithm.fastcore()

    print([i + 1 for i in tissue_reactions])

    names_to_test = [reaction_names[i] for i in fastcore_algorithm.properties['core_idx']]

    with model:
        model.remove_reactions([reaction_names[int(i)] for i in range(s_matrix.shape[1]) if i not in tissue_reactions])
        sol1 = model.optimize()

        print(sol1.to_frame())
        print(sol1.to_frame().loc[pd.Series(names_to_test), :])
