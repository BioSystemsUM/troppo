from cobamp.core.models import ConstraintBasedModel
from random import randint

class Tasks(object):
	def __init__(self, S, lb, ub, rx_names, met_names):
		self.model = ConstraintBasedModel(S, list(*zip(lb, ub)), reaction_names=rx_names, metabolite_names=met_names)


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
			self.task_name = 'task_' + str(randint(1,(2**32)-1))
		else:
			self.task_name = task_name

	
