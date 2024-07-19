tINIT tutorial
==================================================

Imports and Setup
--------------------------------------------------

.. code-block:: python

    import pandas as pd
    import cobra
    import re

    from troppo.omics.readers.generic import TabularReader
    from troppo.methods_wrappers import ModelBasedWrapper
    from troppo.omics.integration import ContinuousScoreIntegrationStrategy
    from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties
..

Initial Setup
--------------------------------------------------

First, we need to define the parsing rules for the GPRs that will be used later on.

.. code-block:: python

    patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
    replace_alt_transcripts = lambda x: patt.sub('', x)
..

Secondly, we need to load the model and the csv file containing the expression dataset.

.. code-block:: python

    model = cobra.io.read_sbml_model('data/HumanGEM_Consistent.xml')
    expression_data = pd.read_csv('data/expression_data.csv', index_col=0)
..

Create the OmicsContainer object
--------------------------------------------------

The `TabularReader` class is used to read and store the omics data in a container that can then be used by *Troppo*.

Relevant arguments from the `TabularReader` class:

- `path_or_df`: the omics data can be either a pandas dataframe or a path to a dataset file. The file can be in any format supported by pandas.
- `index_col`: the name of the column that contains the identifiers of the genes.
- `sample_in_rows`: a boolean indicating whether the samples are in rows or columns.
- `header_offset`: the number of rows to skip before reading the header.
- `omics_type`: a string containing the type of omics data. This is used to select the appropriate integration method.
- `nomenclature`: a string containing the nomenclature of the identifiers in the omics data. This is used to map the identifiers to the identifiers in the model.

The `to_containers()` method returns a list of containers, one for each sample of the dataset. In this example, we will be using only one sample, however, the process can be iterated for all the samples in the dataset.
The `get_integrated_data_map()` method is used to map the identifiers in the omics data to the identifiers in the model. This is done by using the `gpr_gene_parse_function` argument from the `ModelBasedWrapper` class.

.. code-block:: python

    omics_container = TabularReader(path_or_df=omics_data, nomenclature='entrez_id',
                                    omics_type='transcriptomics').to_containers()
    single_sample = omics_container[0]
..

Create a model wrapper
--------------------------------------------------

The `ModelBasedWrapper` class is used to wrap the model so that it can be used by *Troppo*.

Relevant arguments from this class include:

- `model`: the model to be wrapped.
- `ttg_ratio`: the ratio between the number of reactions to be selected and the total number of reactions in the model.
- `gpr_gene_parse_function`: a function that parses the GPRs of the model. This is used to map the identifiers in the omics data to the identifiers in the model.

Important attributes from this class include:

- `model_reader`: a COBRAModelObjectReader instance containing all the information of the model, such as, reaction_ids, metabolite_ids, GPRs, bounds, etc.
- `S`: the stoichiometric matrix of the model.
- `lb`: the lower bounds of the reactions in the model.
- `ub`: the upper bounds of the reactions in the model.

In this specific example we will use the `ReconstructionWrapper` class instead of the base `ModelBasedWrapper` class.

.. code-block:: python

    model_wrapper = ReconstructionWrapper(model=model, ttg_ratio=9999,
                                          gpr_gene_parse_function=replace_alt_transcripts)
..

Map gene IDs in the data to model IDs
---------------------------------------------------

For this we can use the `get_integrated_data_map()` method from the `TabularReader` class. This maps the gene ids in the omics dataset reaction ids in the model through their GPRs, and attributes a score to each reaction in accordance with the expression values of the associated genes. This method returns a dictionary with the reaction ids as keys and the scores as values.

Important arguments from this method include:

- `model_reader`: a COBRAModelObjectReader instance containing all the information of the model. It can be accessed through the `model_wrapper.model_reader`.
- `and_func`: a function that is used to combine the scores of the genes associated with a reaction for AND rules in the GPR. In this example, we will be using the minimum function, which means that the score of a reaction with AND in their GPRs will be the minimum score of the genes associated with it.
- `or_func`: a function that is used to combine the scores of the genes associated with a reaction for OR rules in the GPR. In this example, we will be using the sum function, which means that the score of a reaction with OR in their GPRs will be the sum of the scores of the genes associated with it.

.. code-block:: python

    data_map = single_sample.get_integrated_data_map(model_reader=model_wrapper.model_reader,
                                                     and_func=min, or_func=sum)
..

Integrate Scores
--------------------------------------------------

The `integrate()` method from the `ContinuousScoreIntegrationStrategy` class is used to integrate the scores of the reactions in the model. This method returns a dictionary with the reaction ids as keys and the integrated scores as values. In the case of this continuous scoring method, the resulting scores are the same as the scores in the data map. However, for other scoring methods, such as threshold scoring methods, the result will be a list of reactions with a score above the selected threshold.

Moreover, this method allows us to apply an additional function to the method, which can be useful if you have any protected reactions that need to be in the final model or to remove nan values from the result. This can be done by passing the function as the `score_apply` argument of the `ContinuousScoreIntegrationStrategy` class.

In this example, we will be using a function that replaces the nan values with 0 and returns a list with all the scores. This is the required format for the *tINIT* method.
Keep in mind that if you want to alter this function, the ouput must keep the same format.

.. code-block:: python

    def score_apply(reaction_map_scores):
        return {k:0  if v is None else v for k, v in reaction_map_scores.items()}

    continuous_integration = ContinuousScoreIntegrationStrategy(score_apply=score_apply)
    scores = continuous_integration.integrate(data_map=data_map)
..

For the `ContinuousScoreIntegrationStrategy` the output will be a dictionary with reaction names as keys and scores as values, which is a requirement for the *tINIT* method.
Keep in mind that different integration strategies will have different outputs.
For instance, the `ThresholdScoreIntegrationStrategy` will return a list of reactions with a score above the selected threshold.
Hence, if you use this strategy, *tINIT* will not be able to run.

Run tINIT
--------------------------------------------------

The `tINITProperties` class is used to create the properties for the `tINIT` algorithm. This class contains the following arguments:

- `reactions_scores`: a list of scores for each reaction in the model. This can be obtained from the `integrate()` method of the `ContinuousScoreIntegrationStrategy` class.
- `present_metabolites`: The metabolites that are present in the model.
- `essential_reactions`: The reactions that are essential in the model.
- `allow_excretion`: Whether to allow excretion.
- `no_reverse_loops`: Whether to allow reverse loops.
- `solver`: The solver to be used.

The `tINIT` class is used to run the `tINIT` algorithm. This class contains the following arguments:

- `S`: the stoichiometric matrix of the model. It can be accessed through the `model_wrapper.S`.
- `lb`: the lower bounds of the reactions in the model. It can be accessed through the `model_wrapper.lb`.
- `ub`: the upper bounds of the reactions in the model. It can be accessed through the `model_wrapper.ub`.
- `properties`: a `tINITProperties` instance containing the properties for the `tINIT` algorithm.

In the end, the `run()` method of the `tINIT` class will return a list with indexes of the reactions to be kept in the final model.

.. code-block:: python

    # Create the properties for the GIMME algorithm.
    properties = tINITProperties(reactions_scores=[v for k, v in scores.items()], solver='CPLEX')

    # Run the GIMME algorithm.
    tinit = tINIT(S=model_wrapper.S, lb=model_wrapper.lb, ub=model_wrapper.ub, properties=properties)

    model_tinit = tinit.run()
..
