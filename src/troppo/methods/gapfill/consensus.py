

## TODO: Abstract class for this
class CombinatorialGapfill(object):
	def __init__(self, template_model, task_dict):
		'''
		:param template_model: dict - must contain {S, lb, ub}
		:param task_dict: dict - ??
		'''

		self.template = template_model
		self.tasks = task_dict

	def generate_consensus_model(self, partial_models):
		pass