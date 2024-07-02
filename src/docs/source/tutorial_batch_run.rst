Batch integration of Omics Data
=======================================

Several integration algorithms were introduced in the previous tutorials.
However, the demonstrated approach was limited to a single sample.
In some cases, multiple samples are available and the context-specific models are required for each.
Hence, making the integration of multiple samples a necessity.

`batch_run` is a function from *Cobamp* that allows multiprocessing and is fully compatible with the *Troppo* framework.
Thus allowing the integration of multiple samples in a single run.
This function requires four parameters:

- `function`: the function that will run the reconstruction that needs to be parallelized.
- `sequence`: a list with the containers for each sample.
- `paramargs`: a dictionary with the parameters for the function.
- `threads`: the number of parallel processes to run.

Initial setup
-------------

.. code-block:: python

    import pandas as pd
    import cobra
    import re

    from troppo.omics.readers.generic import TabularReader
    from troppo.methods_wrappers import ReconstructionWrapper
    from cobamp.utilities.parallel import batch_run

    patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
    replace_alt_transcripts = lambda x: patt.sub('', x)

    # load the model
    model = cobra.io.read_sbml_model('data\HumanGEM_Consistent_COVID19_HAM.xml')

    # Create the reconstruction wrapper
    model_wrapper = ReconstructionWrapper(model=model, ttg_ratio=9999,
                                          gpr_gene_parse_function=replace_alt_transcripts)

    # load the data
    omics_data = pd.read_csv(filepath_or_buffer=r'data\Desai-GTEx_ensembl.csv', index_col=0)
    omics_data = omics_data.loc[['Lung_Healthy','Lung_COVID19']]

    # creat omics container
    omics_container = TabularReader(path_or_df=omics_data, nomenclature='entrez_id',
                                    omics_type='transcriptomics').to_containers()
..

Define the function to be parallelized
--------------------------------------

This function uses the `run_from_omics` method from the `ReconstructionWrapper` class. This requires the following parameters:

- `omics_data`: the omics data container for the sample.
- `algorithm`: a string containing the algorithm to use for the reconstruction.
- `and_or_funcs`: a tuple with the functions to use for the AND and OR operations of the GPR.
- `integration_strategy`: a tuple with the integration strategy and the function to apply to the scores.
- `solver`: the solver to use for the optimization.
- `**kwargs`: additional parameters for the reconstruction that are specific to used algorithm.

.. code-block:: python

    def reconstruction_function_gimme(omics_container, parameters: dict):

        def score_apply(reaction_map_scores):
            return {k:0  if v is None else v for k, v in reaction_map_scores.items()}

        flux_threshold, obj_frac, rec_wrapper, method = [parameters[parameter] for parameter in
                                          ['flux_threshold', 'obj_frac', 'reconstruction_wrapper',
                                           'algorithm']]

        reac_ids = rec_wrapper.model_reader.r_ids
        metab_ids = rec_wrapper.model_reader.m_ids
        AND_OR_FUNCS = (min, sum)

        return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method,
                                          and_or_funcs=AND_OR_FUNCS,
                                          integration_strategy=('continuous', score_apply),
                                          solver='CPLEX', obj_frac=obj_frac,
                                          objectives=[{'biomass_human': 1}], preprocess=True,
                                          flux_threshold=flux_threshold, reaction_ids=reac_ids,
                                          metabolite_ids=metab_ids)
..

Considering the function above, the parameters for the reconstruction are defined in a dictionary as follows:

.. code-block:: python

    parameters = {'flux_threshold': 0.8, 'obj_frac': 0.8, 'reconstruction_wrapper': model_wrapper,
                  'algorithm': 'gimme'}
..

Run the batch integration
-------------------------

.. code-block:: python

    batch_gimme_res = batch_run(reconstruction_function_gimme, omics_container, parameters, threads=2)

..
