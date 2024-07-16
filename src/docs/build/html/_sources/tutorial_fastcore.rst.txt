FastCORE tutorial
==================================================

Imports and Setup
--------------------------------------------------

.. code-block:: python

    import pandas as pd
    import cobra
    import re

    from troppo.omics.readers.generic import TabularReader
    from troppo.methods_wrappers import ModelBasedWrapper
    from troppo.omics.integration import CustomSelectionIntegrationStrategy
    from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties
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

The `CustomSelectionIntegrationStrategy` class allows the user to define a custom function that is tailored to the output we need for the following steps of the pipeline.
For the `FastCORE` method the most adequate integration strategy is to select reactions whose score is above a defined threshold.

This also can be achieved by using the `ThresholdSelectionIntegrationStrategy` class, however, since we also want to include a set of reaction to be protected during the integration we will use a custom method that will be defined by `integration_fx`.

Moreover, through this function we also want the output to be a list with reaction IDs that will belong to the core reactions that will be inputted for the `FastCORE` algorithm.

.. code-block:: python

    from math import log
    threshold =  (5 * log(2))
    protected_reactions = ['biomass']

    def integration_fx(reaction_map_scores):
        return [[k for k, v in reaction_map_scores.get_scores().items() if (v is not None and v > threshold) or k in protected_reactions]]

    threshold_integration = CustomSelectionIntegrationStrategy(group_functions=[integration_fx])
    threshold_scores = threshold_integration.integrate(data_map=data_map)

    print(threshold_scores)
..

Run FASTcore
--------------------------------------------------

The `FastcoreProperties` class is used to create the properties for the GIMME algorithm. This class contains the following arguments:

- `core`: List of indexes of the reactions that are considered core, as determined by the integrated scores.
- `flux_threshold`: Flux threshold for the algorithm.
- `solver`: Solver to be used.

The `FASTcore` class is used to run the GIMME algorithm. This class contains the following arguments:

- `S`: the stoichiometric matrix of the model. It can be accessed through the `model_wrapper.S`.
- `lb`: the lower bounds of the reactions in the model. It can be accessed through the `model_wrapper.lb`.
- `ub`: the upper bounds of the reactions in the model. It can be accessed through the `model_wrapper.ub`.
- `properties`: a `FastcoreProperties` instance containing the properties for the GIMME algorithm.

In the end, the `run()` method of the `FASTcore` class will return a list with the index of the reactions to be keept in the model.

.. code-block:: python

    # Get the index of the reaction of the CORE reaction set
    ordered_ids = {r:i for i,r in enumerate(model_wrapper.model_reader.r_ids)}
    core_idx = [[ordered_ids[k] for k in l] for l in threshold_scores]

    # Define the FastCORE properties
    properties = FastcoreProperties(core=core_idx, solver='CPLEX')

    # instantiate the FastCORE class
    fastcore = FASTcore(S=model_wrapper.S, lb=model_wrapper.lb, ub=model_wrapper.ub, properties=properties)

    # Run the algorithm
    model_fastcore = fastcore.run()
..
