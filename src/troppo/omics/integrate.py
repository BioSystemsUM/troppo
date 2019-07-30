"""
 Created by Jorge Gomes on 05/04/2018
 TsmRec
 integrate
 
"""
from .omics_data_map import OmicsDataMap

"""
 Function responsible for the integration of different omics data with a metabolic model loaded with framed package.
 see help(integrateOmics)
"""

# Nested function implementation to ease variable scoping


def integrateOmics(model_reader, omics_container, and_func=min, or_func=max, apply_fx=None):
    """
    Function responsible for the integration of different omics data with a metabolic model loaded with framed package.
    Matches model ids for gene_ids, metabolites or reaction ids with those present in the omicsContainer object.

    :param CBModel: (obj) a metabolic model object previously loaded with framed package
    :param omicsContainer: (obj) an omics container object previously created using OmicsContainer class.
    :param modelField: (str) the model field where ids shall be retrieved and used for the integration.(must be either
                        "id" or "name")
    :param and_func:(func) the mathematical function to replace the "AND" operator present in the Gene-Protein-Rules
    :param or_func:(func) the mathematical function to replace the "OR" operator present in the Gene-Protein-Rules


    :return m: (obj) an OmicsDataMap object which contains the mapping between reactions/metabolites and its fluxes
    based on the supplied omics data.
   """

    def g2rIntegrate():
        """
        Handles integration of both proteomics and transcriptomics expression data relying on framed's gene2reaction
        """
        #suffixAndPrefix()
        d = model_reader.g2rx(omics_container.get_Data(), or_fx=or_func, and_fx=and_func, apply_fx=apply_fx)
        return aux_createMap(d, 'ReactionDataMap')

    def aux_createMap(mMap, mapType):
        m = OmicsDataMap(mMap, mapType)
        return m

    # execution commands

    omicsType = omics_container.otype.lower()

    if omicsType.lower() in ['proteomics', 'transcriptomics']:
        return g2rIntegrate()
    else:
        raise Exception('Omics data type not yet supported')








