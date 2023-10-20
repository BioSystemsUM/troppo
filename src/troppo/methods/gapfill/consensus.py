from typing import Iterable

from numpy import array, where, logical_xor
from troppo.tasks.core import TaskEvaluator
from time import time
import warnings


## TODO: Abstract class for this
class CombinatorialGapfill(object):
    """
    This Class uses a combinatorial gap-filling algorithm to remove reactions to the final model

    Parameters
    ----------
    template_model : dict
        A dictionary containing the template model. Must contain {S, lb, ub} and possibly rx_names/met_names.
    task_dict : dict
        A dictionary containing the tasks that need to be completed by the final model.
    min_acceptable_tasks : int or float (default: 0.5)
        The minimum number of tasks that need to be completed by the final model. If an integer is provided,
        this is the minimum number of tasks that need to be completed. If a float is provided, this is the
        minimum percentage of tasks.
    """

    def __init__(self, template_model, tasks, min_acceptable_tasks=0.5):
        self.template = template_model
        self.tasks = tasks

        if isinstance(min_acceptable_tasks, int):
            self.__min_tasks = min_acceptable_tasks

        elif isinstance(min_acceptable_tasks, float) and (0 < min_acceptable_tasks <= 1):
            self.__min_tasks = (min_acceptable_tasks // len(self.tasks)) if (len(self.tasks) > 0) else 0

        else:
            raise ValueError('Invalid parameter `min_acceptable_tasks`')

    @staticmethod
    def generate_partial_models(rx_presences: list) -> dict:
        """
        Generates partial models from the presence/absence of reactions in the template model.

        Parameters
        ----------
        rx_presences: list
            A list containing the indices of the reactions that are present in the model.

        Returns
        -------
        partial_models: dict
            A dictionary containing the partial models.
        """
        rp = array(rx_presences).astype(int).sum(axis=0)

        partial_models = {i + 1: where(rp > i + 1)[0] for i in range(len(rx_presences))}
        partial_models[0] = where(array([True] * len(rp)))[0]

        return partial_models

    def __get_task_validator_instances(self, partial_models: dict):
        """
        Creates a list of TaskValidator instances for each partial model.

        Parameters
        ----------
        partial_models: dict
            A dictionary containing the partial models.

        Returns
        -------
        task_instances: dict
            A list of TaskValidator instances.
        """
        task_instances = {}

        for i, items in partial_models.items():
            model_params = {k: v for k, v in self.template.items()}
            model_params['context'] = [model_params['rx_names'][k] for k in
                                       items] if 'rx_names' in model_params else items

            task_instances[i] = TaskEvaluator(**model_params)

        return task_instances

    def gapfill(self, rx_presences: list) -> Iterable:
        """
        This method performs gap-filling to build the final model.

        Parameters
        ----------
        rx_presences: list
            A list containing the indices of the reactions that are present in the model.

        Returns
        -------
        final_model: list
            A list containing the indices of the reactions that are present in the final model.

        """
        print('Generating partial models...')
        partials = self.generate_partial_models(rx_presences)

        print('Creating validator instances...')
        validators = self.__get_task_validator_instances(partials)

        satisfied, j = True, 0
        completed_tasks, lost_metabolic_tasks, lost_reaction_set = {}, {}, {}

        print('Validating partial model tasks...')
        while (j < len(validators)) and satisfied:
            print('\tModel', j, 'with', len(partials[j]), 'reactions')
            valid_tasks = set([k for k, v in validators[j].evaluate(self.tasks).items() if v])
            completed_tasks[j] = valid_tasks

            if (j > 0) and (j < len(rx_presences)):
                lost_metabolic_tasks[j - 1] = completed_tasks[j - 1] - completed_tasks[j]
                lost_reaction_set[j - 1] = logical_xor(rx_presences[j], rx_presences[j - 1])
                lost_reaction_set[j - 1] = set(where(lost_reaction_set[j - 1])[0])

            satisfied = len(valid_tasks) >= self.__min_tasks
            print('\tModel', j, 'completes', len(valid_tasks), 'out of ', len(self.tasks))

            if satisfied:
                j += 1

        # Step 4
        validators = {i: validators[i] for i in range(j)}

        print('Building final model...')
        return self.build_final_model(partials={i: partials[i] for i in range(j)},
                                      start_from=j - 1,
                                      lost_metabolic_tasks=lost_metabolic_tasks,
                                      lost_reaction_sets=lost_reaction_set,
                                      validators=validators)

    def build_final_model(self, partials: dict, start_from: int, lost_metabolic_tasks: dict,
                          lost_reaction_sets: dict, validators: dict) -> set:
        """
        This method builds the final model by adding reactions to the partial model that completes the most tasks.

        Parameters
        ----------
        partials: dict
            A dictionary containing the partial models.
        start_from: int
            The index of the partial model to start from.
        lost_metabolic_tasks: dict
            A dictionary containing the tasks that were lost when a reaction was removed from the model.
        lost_reaction_sets: dict
            A dictionary containing the reactions that were removed from the model.
        validators: dict
            A dictionary containing the TaskValidator instances.

        Returns
        -------
        final_model: set
            A set containing the indices of the reactions that are present in the final model.

        """
        to_del = set()
        final_model = set()

        for i in range(start_from, -1, -1):
            final_model, model_validator = set(partials[i]), validators[i]
            tasks, reactions = lost_metabolic_tasks[i], set(lost_reaction_sets[i])
            test_kos = to_del | set(reactions)

            print('\tTesting model', i, 'with', len(test_kos), 'reaction knockouts.')

            times_spent = []

            if len(tasks) >= 1:
                for r in test_kos:
                    start_time = time()
                    valid = self.is_valid_model(tasks, r, model_validator)
                    if valid:
                        final_model -= {r}
                        to_del |= {r}
                        if i > 0:
                            lost_reaction_sets[i - 1] |= {r}

                    end_time = time()
                    times_spent.append(end_time - start_time)

            else:
                final_model -= test_kos
                to_del |= test_kos
                warnings.warn(' '.join(['Model', str(i), 'not tested.', 'No lost metabolic tasks']))

            if len(times_spent) > 0:
                print('\t', (sum(times_spent) / len(times_spent)) * 1000, 'ms per optimization')

        return final_model

    def is_valid_model(self, tasks: set, ko: int, validator: TaskEvaluator) -> bool:
        """
        This method checks if a model is valid by testing the presence of the tasks after removing a reaction.

        Parameters
        ----------
        tasks: set
            A set containing the indices of the tasks to be tested.
        ko: int
            The index of the reaction to be knocked out.
        validator: TaskEvaluator
            A TaskEvaluator instance.

        Returns
        -------
        valid: bool
            A boolean indicating if the model is valid.
        """
        lb, ub = validator.model.get_reaction_bounds(ko)
        validator.model.set_reaction_bounds(ko, lb=0, ub=0)

        task_validation = validator.evaluate([k for k in self.tasks if k.task_name in tasks])
        validator.model.set_reaction_bounds(ko, lb=lb, ub=ub)

        return not (False in task_validation.values())
