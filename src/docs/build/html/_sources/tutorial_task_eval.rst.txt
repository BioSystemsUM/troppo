Task performance evaluation
=======================================================================

*Troppo* includes a module that allows for the evaluation of metabolic tasks performance.
The class responsible for this is called `TaskEvaluator`.
This class is a wrapper around a model that allows the evaluation of tasks on the model. It can be used to evaluate a single task, or to evaluate a batch of tasks on a batch of models.

To initialize a `TaskEvaluator` object, you need to pass a model object and the tasks to evaluate.
The tasks should be instances of the `Task` class, which is a simple data structure that contains the following fields:

- `reaction_dict`: a dictionary with the reaction identifiers as keys and a dictionary with the metabolites and their respective stoichiometry as values. (eg. rxd = {'r1':({'m1':-1, 'm2':2}, (lb, ub)), ... })
- `inflow_dict`: a dictionary with the metabolite identifiers as keys and the inflow rate as values. (eg. inflow = {'m1':(1, 1), ... })
- `outflow_dict`: a dictionary with the metabolite identifiers as keys and the outflow rate as values. (eg. outflow = {'m5':(5, 5), ... })

Setup the model and the tasks
_______________________________________________________________________

.. code-block:: python

    #imports
    from troppo.tasks.core import TaskEvaluator
    from troppo.tasks.task_io import JSONTaskIO
    import pandas as pd
    from json import JSONEncoder, JSONDecoder
    from cobamp.utilities.parallel import batch_run

    from numpy import log
    import re

    #load the model
    task_model = read_sbml_model('data/Recon3D_301_consistent.xml')

    #load the reconstruction results as a dictionary
    fastcore_res_dict = pd.read_csv('data/r3d_compact_ccle_bc_fastcore.csv',
                                    index_col=[0, 1]).T.to_dict()

    #get only first sample
    sample = list(fastcore_res_dict.keys())[0]
    fastcore_res_dict = {sample: fastcore_res_dict[sample]}

    TASKS_PATH = 'data/nl2019_tasks_r3d_compact.json'

    # parse tasks from a previously existing JSON
    # the supplied file contains tasks adapted from the publication of Richelle et. al, 2019
    task_list = [t for t in JSONTaskIO().read_task(TASKS_PATH) if len((set(t.inflow_dict) |
                 set(t.outflow_dict)) - set([m.id for m in task_model.metabolites])) == 0]

    for task in task_list:
        task.inflow_dict = {k: v if k not in task.outflow_dict.keys() else [-1000, 1000] for k, v in
                            task.inflow_dict.items()}
        task.outflow_dict = {k: v for k, v in task.outflow_dict.items()
                             if k not in task.inflow_dict.items()}
    for task in task_list:
        task.mandatory_activity = []

    # tasks should be evaluated without open boundary reactions. We can easily close them on
    # the COBRA model
    for k in task_model.boundary:
        k.knock_out()

    # get the names of all reactions in the model - this will be useful further on
    all_reactions = set([r.id for r in task_model.reactions])
..

Evaluate a set of tasks
_______________________________________________________________________

.. code-block:: python

    # create a structure to hold all of these results - a dictionary
    task_eval_results = {}

    # for each k (tuple with algorithm and sample information) and result (dict with reaction presences)
    for k, result in fastcore_res_dict.items():
        # using with statements to change the COBRA model temporarily
        # this is done to knock-out reaction not appearing the FASTCORE result
        with task_model as context_specific_model:
            # get reactions included in the sample-specific model
            protected = set([k for k, v in result.items() if v])
            # get reactions except the protected ones
            to_remove = all_reactions - protected
            for rid in to_remove:
                # knock-out reactions not in the model
                context_specific_model.reactions.get_by_id(rid).knock_out()

            # create a task evaluator instance with the context specific model and the supplied
            # task list and solver
            task_eval = TaskEvaluator(model=context_specific_model, tasks=task_list, solver='CPLEX')

            # get task names (for future reference)
            task_names = task_eval.tasks

            # use the batch_function from the TaskEvaluator class (takes the name of a loaded task,
            # a params dictionary with the task evaluator associated to the 'tev' key) and set the
            # amount of threads to be used
            batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names,
                                        {'tev': task_eval}, threads=1)
        # each element in the list of results in batch_res_tasks is a tuple of length 3:
        # 0 - boolean flag representing the task evaluation
        # 1 - Solution instance used to evaluate the task
        # 2 - A dictionary with reactions supposed to be active mapped to True/False
        # according to that criterion

        # keep only items 0 and 2 of the task result - we don't need the flux distribution
        task_csm_res = {k: (v[0], v[2]) for k, v in dict(zip(task_names, batch_res_tasks)).items()}
        print(k, len(protected), len([v for k, v in task_csm_res.items() if v[0]]), 'tasks completed.')
        # assign this dictionary to it's sample on the master results dictionary
        task_eval_results[k] = task_csm_res
..

Save the results
_______________________________________________________________________

.. code-block:: python

    # save these results for later analysis as a JSON file
    with open('data/r3d_compact_task_results_ccle_bc_new_nodrains_only_feas.json', 'w') as f:
        f.write(JSONEncoder().encode([(k, v) for k, v in task_eval_results.items()]))