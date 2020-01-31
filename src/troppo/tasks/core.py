from random import randint

from pathos.multiprocessing import cpu_count
from cobamp.utilities.parallel import batch_run
from cobamp.utilities.context import CommandHistory
from cobamp.core.models import ConstraintBasedModel
from cobamp.core.optimization import Solution
from cobamp.wrappers.external_wrappers import get_model_reader

MP_THREADS = cpu_count()
INF = float('inf')

from functools import partial

#
# def disable_task_bounds(model, task_reactions):
# 	for rx in task_reactions:
# 		model.set_reaction_bounds(rx, lb=0, ub=0)
#
# def enable_task_bounds(model, task_items):
# 	## TODO: excepção
# 	for d in task_items:
# 		for rx, bounds in d.items():
# 			lb, ub = bounds
# 			model.set_reaction_bounds(rx, lb=lb, ub=ub)
#
# def _init_task_solver(model, task_components, task_fail_status, flux_constraints, task_added_reactions):
# 	global _model, _task_components, _task_fail_status, _flux_constraints, _task_added_reactions
# 	_model, _task_components, _task_fail_status, _flux_constraints, _task_added_reactions = \
# 		model, task_components, task_fail_status, flux_constraints, task_added_reactions
#
# def _task_iteration(task_name):
# 	global _model, _task_components, _task_fail_status, _flux_constraints, _task_added_reactions
# 	return task_name, evaluate_task(_model, _task_components[task_name], _task_fail_status[task_name],
# 									_flux_constraints[task_name], _task_added_reactions)
#
# def evaluate_task(model, task_items, sfail, flux_constraints, task_reactions):
# 	enable_task_bounds(model, task_items)
# 	sol = model.optimize()
# 	constraint_compliant = evaluate_flux_constraints(flux_constraints, sol.var_values())
# 	is_optimal, is_infeasible = [sol.status() == x for x in ['optimal', 'infeasible']]
# 	task_status = (not is_optimal and is_infeasible) if sfail else (is_optimal and not is_infeasible)
# 	if task_status:
# 		task_status = (task_status or (not constraint_compliant)) if sfail else (task_status and constraint_compliant)
# 	disable_task_bounds(model, task_reactions)
# 	return task_status
#
# def evaluate_flux_constraints(flux_constraints, solution):
# 	for i, v in flux_constraints.items():
# 		lb,ub = [(INF*sign) if (k == None) else k for k,sign in zip(v,[-1, 1])]
# 		valid = lb <= solution[i] <= ub
# 		if not valid:
# 			return False
# 	return True
#
# def task_pool(model, task_components, task_fail_status, flux_constraints, task_added_reactions):
# 	threads = MP_THREADS
# 	task_list = list(task_components.keys())
# 	res_map = {r: i for i, r in enumerate(task_list)}
# 	true_threads = min((len(task_list) // 2) + 1, threads)
# 	it_per_job = len(task_list) // threads
# 	pool = _ProcessPool(
# 		processes=true_threads,
# 		initializer=_init_task_solver,
# 		initargs=(model, task_components, task_fail_status, flux_constraints, task_added_reactions)
# 	)
# 	for i, value in pool.imap_unordered(_task_iteration, task_list,
# 										chunksize=it_per_job):
# 		res_map[i] = value
#
# 	pool.close()
# 	pool.join()
#
# 	return res_map
#
# ## TODO: For now, context will be applied using knockouts. In the future we should actually prune and simplify
#
# def apply_context(model, context):
# 	for i,index in enumerate(model.reaction_names):
# 		if (index not in context) and (i not in context):
# 			model.set_reaction_bounds(index, lb=0, ub=0)
#
#
# class TaskEvaluator(object):
# 	def __init__(self, S, lb, ub, rx_names, met_names, objective_reaction=None, context=None):
# 		self.model = ConstraintBasedModel(S, list(zip(lb, ub)), reaction_names=rx_names, metabolite_names=met_names)
# 		if context != None:
# 			apply_context(self.model, context)
# 		self.tasks = {}
# 		self.task_fail_status = {}
# 		self.__add_reactions = []
# 		self.flux_constraints = {}
#
# 		if objective_reaction:
# 			self.model.set_objective({objective_reaction: 1}, False)
# 		else:
# 			self.model.set_objective({0:1}, False)
#
# 	def populate_model(self, task):
# 		rd, ifd, ofd, flc = task.get_task_components()
# 		sfail = task.should_fail
#
# 		## TODO: excepções para qd um metabolito não existir
# 		task_equations = {}
# 		for r_id, tup in rd.items():
# 			r_mets, bounds = tup
# 			r_name = '_'.join([task.task_name, r_id])
# 			self.model.add_reaction(arg=r_mets, bounds=bounds, name=r_name)
# 			task_equations[r_name] = bounds
# 			self.__add_reactions.append(r_name)
#
# 		task_sinks = {}
# 		for coef, sinkdict in zip([-1, 1], [ofd, ifd]):
# 			for m_id, bounds in sinkdict.items():
# 				r_name = '_'.join([task.task_name, m_id, 'sink'])
# 				self.model.add_reaction(arg={m_id:coef}, bounds=bounds, name=r_name)
# 				self.__add_reactions.append(r_name)
# 				task_sinks[r_name] = bounds
# 		self.tasks[task.task_name] = (task_equations, task_sinks)
# 		self.task_fail_status[task.task_name] = sfail
# 		self.flux_constraints[task.task_name] = flc
#
# 	def __run_single_task(self, task):
# 		return evaluate_task(model=self.model,
# 							 task_items=self.tasks[task],
# 							 sfail=self.task_fail_status[task],
# 							 flux_constraints=self.flux_constraints[task],
# 							 task_reactions=self.__add_reactions)
#
# 	def __task_initializer(self, task):
# 		if task.task_name not in self.tasks.keys():
# 			self.populate_model(task)
#
# 	def evaluate(self, task_arg):
# 		if isinstance(task_arg, (list, tuple, set)):
# 			for task in task_arg:
# 				if isinstance(task, Task):
# 					self.__task_initializer(task)
# 				else:
# 					raise TypeError('Invalid type object found within the task_arg iterable. Expected Task, found'+str(type(task)))
# 			if len(task_arg) <= MP_THREADS:
# 				return {task.task_name:self.__run_single_task(task.task_name) for task in task_arg}
# 			else:
# 				return task_pool(self.model, self.tasks, self.task_fail_status, self.flux_constraints, self.__add_reactions)
#
# 		elif isinstance(task_arg, Task):
# 			self.__task_initializer(task_arg)
# 			return self.__run_single_task(task_arg.task_name)
# 		else:
# 			raise TypeError('task_arg expected to be of type Task. Found '+str(type(task_arg))+' instead')

class TaskEvaluator(object):
	def __init__(self, **kwargs):
		if 'solver' in kwargs:
			solver = kwargs['solver']
		else:
			solver = None

		if 'model' in kwargs.keys():
			model_obj = kwargs['model']
			if isinstance(model_obj, ConstraintBasedModel):
				self.model = model_obj
			else:
				self.model = get_model_reader(model_obj).to_cobamp_cbm(solver)
		else:
			if 'lb' in kwargs.keys():
				S, lb, ub, rxn, mtn = [kwargs[k] for k in ['S','lb','ub','reaction_names','metabolite_names']]
				bounds = list(zip(lb, ub))

				self.model = ConstraintBasedModel(S, bounds, rxn, mtn, True, solver)

	def evaluate_task(self, task):
		pass

class Task(object):

	__defaults__ = {
		'should_fail': False,
		'reaction_dict': {},
		'flow_dict': {},
		'flux_constraints': {},
		'mandatory_activity': [],
		'name': 'default_task',
		'annotations': {}
	}

	__types__ = {
		'reaction_dict': dict,
		'flow_dict': dict,
		'should_fail': bool,
		'name': str,
		'flux_constraints': dict,
		'annotations': dict,
		'mandatory_activity': list
	}

	def __init__(self, **kwargs):

		'''
		reaction_dict: rxd = {'r1':({'m1':-1, 'm2':2}, (lb, ub)), ...
		inflow_dict: ifd = {'m3':(1,1), ...
		outflow_dict: ofd = {'m5':(5,5), ...

		'''
		for k,v in self.__defaults__.items():
			itype, dval = self.__types__[k], self.__defaults__[k]
			setattr(self, k, dval if k not in kwargs.keys() else kwargs[k])

	@property
	def should_fail(self):
		return self.__should_fail

	@should_fail.setter
	def should_fail(self, value):
		self.__should_fail = value

	@property
	def reaction_dict(self):
		return self.__reaction_dict

	@reaction_dict.setter
	def reaction_dict(self, value):
		self.__reaction_dict = value

	@property
	def flux_constraints(self):
		return self.__flux_constraints

	@flux_constraints.setter
	def flux_constraints(self, value):
		self.__flux_constraints = value

	@property
	def flow_dict(self):
		return self.__flow_dict

	@flow_dict.setter
	def flow_dict(self, value):
		self.__flow_dict = value

	@property
	def name(self):
		return self.__name

	@name.setter
	def name(self, value):
		self.__name = value

	@property
	def annotations(self):
		return self.__annotations

	@annotations.setter
	def annotations(self, value):
		self.__annotations = value

	@property
	def mandatory_activity(self):
		return self.__mandatory_activity

	@mandatory_activity.setter
	def mandatory_activity(self, value):
		self.__mandatory_activity = value

	def get_task_manipulation_cmds(self, model: ConstraintBasedModel):
		## reaction_dict - add reactions to the model

		command_history = CommandHistory()
		for k,v in self.reaction_dict.items():
			reaction_name = '_'.join([self.name,k,'flow'])
			command_history.queue_command(model.add_reaction, {'arg':v[0], 'bounds':v[1], 'name': reaction_name})

		## flow_dict - add drains to the model
		for k,v in self.flow_dict.items():
			sink_name = '_'.join([self.name,k,'flow'])
			command_history.queue_command(model.add_reaction, {'arg':{k: 1}, 'bounds':v, 'name': sink_name})

		## constraint_dict - impose additional bounds
		for k,v in self.flux_constraints.items():
			command_history.queue_command(model.set_reaction_bounds, {'index':k, 'lb':v[0], 'ub':v[1]})

		return command_history

	def apply_evaluate(self, model: ConstraintBasedModel):
		with model as task_model:
			commands = self.get_task_manipulation_cmds(task_model)
			commands.execute_all(True)
			task_evaluation = self.evaluate_solution(task_model.optimize())
		return task_evaluation

	def evaluate_solution(self, sol: Solution, ftol=1e-6):
		is_optimal = sol.status() == 'optimal'
		mandatory_are_valid = sum([abs(sol[k]) > ftol for k in self.mandatory_activity]) == len(self.mandatory_activity)
		conditions_are_met = (is_optimal and mandatory_are_valid)

		return (conditions_are_met and not self.should_fail)



if __name__ == '__main__':
	from numpy import array

	S = array([[1, -1, 0, 0, -1, 0, -1, 0, 0],
			   [0, 1, -1, 0, 0, 0, 0, 0, 0],
			   [0, 1, 0, 1, -1, 0, 0, 0, 0],
			   [0, 0, 0, 0, 0, 1, -1, 0, 0],
			   [0, 0, 0, 0, 0, 0, 1, -1, 0],
			   [0, 0, 0, 0, 1, 0, 0, 1, -1]])

	rx_names = ["R" + str(i) for i in range(1, 10)]
	mt_names = ["M" + str(i) for i in range(1, 7)]

	irrev = [0, 1, 2, 4, 5, 6, 7, 8]
	bounds = [(0 if i in irrev else -1000, 1000) for i in range(9)]
	lb, ub = list(zip(*bounds))
	T = array([0] * S.shape[1]).reshape(1, S.shape[1])
	T[0, 8] = -1
	b = array([-1]).reshape(1, )

	# tasks = TaskEvaluator(S, lb, ub, rx_names, mt_names)

	task1 = Task(
		flow_dict={'M4': (-4, -4)},
		should_fail=True,
		flux_constraints={'R2': (6, 10)})

	from cobamp.core.models import ConstraintBasedModel

	cbm = ConstraintBasedModel(S, list(zip(lb, ub)), rx_names, mt_names)
	task1.apply_evaluate(cbm)
