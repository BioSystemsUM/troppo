from cobamp.core.models import ConstraintBasedModel
from random import randint

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

	def disable_task_bounds(self):
		for rx in self.__add_reactions:
			self.model.set_reaction_bounds(rx, lb=0, ub=0)

	def enable_task_bounds(self, task):
		## TODO: excepção
		for d in self.tasks[task]:
			for rx, bounds in d.items():
				lb, ub = bounds
				self.model.set_reaction_bounds(rx, lb=lb, ub=ub)

	def run_tasks(self, tasks):
		for task in tasks:
			if task not in self.tasks:
				self.populate_model(task)

		self.disable_task_bounds()
		task_sols = {}

		for task_name in self.tasks:
			sfail = self.task_fail_status[task_name]
			self.enable_task_bounds(task_name)
			sol = self.model.optimize()
			is_optimal, is_infeasible = [sol.status() == x for x in ['optimal', 'infeasible']]
			task_status = (not is_optimal and is_infeasible) if sfail else (is_optimal and not is_infeasible)
			task_sols[task_name] = task_status

		return task_sols

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

	task1 = Task(
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
	)

	res = tasks.run_tasks([task1])

	print(res)