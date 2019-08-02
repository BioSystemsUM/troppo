from cobamp.core.models import ConstraintBasedModel
from random import randint
from pathos.pools import _ProcessPool
from pathos.multiprocessing import cpu_count

MP_THREADS = cpu_count()


def disable_task_bounds(model, task_reactions):
	for rx in task_reactions:
		model.set_reaction_bounds(rx, lb=0, ub=0)

def enable_task_bounds(model, task_items):
	## TODO: excepção
	for d in task_items:
		for rx, bounds in d.items():
			lb, ub = bounds
			model.set_reaction_bounds(rx, lb=lb, ub=ub)

def _init_task_solver(model, task_components, task_fail_status, task_added_reactions):
	global _model, _task_components, _task_fail_status, _task_added_reactions
	_model, _task_components, _task_fail_status, _task_added_reactions = model, task_components, task_fail_status, task_added_reactions

def _task_iteration(task_name):
	global _model, _task_components, _task_fail_status, _task_added_reactions
	return task_name, evaluate_task(_model, _task_components[task_name], _task_fail_status[task_name], _task_added_reactions)

def evaluate_task(model, task_items, sfail, task_reactions):
	enable_task_bounds(model, task_items)
	sol = model.optimize()
	is_optimal, is_infeasible = [sol.status() == x for x in ['optimal', 'infeasible']]
	task_status = (not is_optimal and is_infeasible) if sfail else (is_optimal and not is_infeasible)
	disable_task_bounds(model, task_reactions)
	return task_status

def task_pool(model, task_components, task_fail_status, task_added_reactions):
	threads = MP_THREADS
	task_list = list(task_components.keys())
	res_map = {r: i for i, r in enumerate(task_list)}
	true_threads = min((len(task_list) // 2) + 1, threads)
	it_per_job = len(task_list) // threads
	pool = _ProcessPool(
		processes=true_threads,
		initializer=_init_task_solver,
		initargs=(model, task_components, task_fail_status, task_added_reactions)
	)
	for i, value in pool.imap_unordered(_task_iteration, task_list,
										chunksize=it_per_job):
		res_map[i] = value

	pool.close()
	pool.join()

	return res_map



class Tasks(object):
	def __init__(self, S, lb, ub, rx_names, met_names, objective_reaction=None):
		self.model = ConstraintBasedModel(S, list(zip(lb, ub)), reaction_names=rx_names, metabolite_names=met_names)
		self.tasks = {}
		self.task_fail_status = {}
		self.__add_reactions = []

		if objective_reaction:
			self.model.set_objective({objective_reaction: 1}, False)
		else:
			self.model.set_objective({0:1}, False)

	def populate_model(self, task):
		rd, ifd, ofd = task.get_task_components()
		sfail = task.should_fail

		## TODO: excepções para qd um metabolito não existir
		task_equations = {}
		for r_id, tup in rd.items():
			r_mets, bounds = tup
			r_name = '_'.join([task.task_name, r_id])
			self.model.add_reaction(arg=r_mets, bounds=bounds, name=r_name)
			task_equations[r_name] = bounds
			self.__add_reactions.append(r_name)

		task_sinks = {}
		for coef, sinkdict in zip([-1, 1], [ofd, ifd]):
			for m_id, bounds in sinkdict.items():
				r_name = '_'.join([task.task_name, m_id, 'sink'])
				self.model.add_reaction(arg={m_id:coef}, bounds=bounds, name=r_name)
				self.__add_reactions.append(r_name)
				task_sinks[r_name] = bounds
		self.tasks[task.task_name] = (task_equations, task_sinks)
		self.task_fail_status[task.task_name] = sfail

	def __run_single_task(self, task):
		return evaluate_task(self.model, self.tasks[task], self.task_fail_status[task], self.__add_reactions)

	def __task_initializer(self, task):
		if task.task_name not in self.tasks.keys():
			self.populate_model(task)

	def evaluate(self, task_arg):
		if isinstance(task_arg, (list, tuple, set)):
			for task in task_arg:
				if isinstance(task, Task):
					self.__task_initializer(task)
				else:
					raise TypeError('Invalid type object found within the task_arg iterable. Expected Task, found'+str(type(task)))
			if len(task_arg) <= MP_THREADS:
				return {task.task_name:self.__run_single_task(task.task_name) for task in task_arg}
			else:
				return task_pool(self.model, self.tasks, self.task_fail_status, self.__add_reactions)

		elif isinstance(task_arg, Task):
			self.__task_initializer(task_arg)
			return self.__run_single_task(task_arg.task_name)
		else:
			raise TypeError('task_arg expected to be of type Task. Found '+str(type(task_arg))+' instead')





class Task(object):
	def __init__(self, reaction_dict, inflow_dict, outflow_dict, should_fail, task_name=None):
		'''
		reaction_dict: rxd = {'r1':({'m1':-1, 'm2':2}, (lb, ub)), ...
		inflow_dict: ifd = {'m3':(1,1), ...
		outflow_dict: ofd = {'m5':(5,5), ...

		'''
		self.reaction_dict = reaction_dict
		self.inflow_dict = inflow_dict
		self.outflow_dict = outflow_dict
		self.should_fail = should_fail

		if (task_name == None) and (not isinstance(task_name, str)):
			self._task_name = 'task_' + str(randint(1,(2**32)-1))
		else:
			self._task_name = task_name


	def get_task_components(self):
		return self.reaction_dict, self.inflow_dict, self.outflow_dict

	@property
	def task_name(self):
		return self._task_name

	def should_fail(self):
		return self.should_fail



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

	tasks = Tasks(S, lb, ub, rx_names, mt_names)

	task1 = [Task(
		reaction_dict={
			'A': ({'M1':-1, 'M2':2},(10,10)),
			'B': ({'M3':-2, 'M6':3},(2,10))
		},
		inflow_dict={
			'M2':(10000,10000)
		},
		outflow_dict={
			'M4':(4,10)
		},
		should_fail=True
	) for i in range(20)]

	res = tasks.evaluate(task1)

	print(res)