import abc
from io import TextIOBase
from json import JSONDecoder, JSONEncoder
import warnings
from troppo.tasks.core import Task

class TaskIO(object):
	__metaclass__ = abc.ABCMeta

	def read_task(self, buffer_or_path):
		if isinstance(buffer_or_path, TextIOBase):
			return self.read_from_string(buffer_or_path.read())

		elif isinstance(buffer_or_path, str):
			with open(buffer_or_path, 'r') as f:
				return self.read_from_string(f.read())

		else:
			raise TypeError('Invalid buffer or path')

	def write_task(self, buffer_or_path, task_arg):
		task_string = self.write_to_string(task_arg)
		if isinstance(buffer_or_path, TextIOBase):
			buffer_or_path.write(task_string)
		elif isinstance(buffer_or_path, str):
			with open(buffer_or_path, 'w') as f:
				return f.write(task_string)
		else:
			raise TypeError('Invalid buffer or path')

	@abc.abstractmethod
	def write_to_string(self, task):
		return ''

	@abc.abstractmethod
	def read_from_string(self, string):
		return []


class JSONTaskIO(TaskIO):
	def read_from_string(self, string):
		def sanity_check(json_dict):
			types, defaults = Task.__types__, Task.__defaults__
			for key in defaults:
				if key not in json_dict.keys():
					str_msg = ' '.join(['Key',key,'has no value.','Setting default value =',str(defaults[key])])
					warnings.warn(str_msg)
					json_dict[key] = defaults[key]

				elif type(json_dict[key]) != types[key]:
					str_msg = ' '.join(['Key',key,'with value=',str(json_dict[key]),'does not match the expected type.','Setting default value =',str(defaults[key])])
					warnings.warn(str_msg)
					json_dict[key] = defaults[key]
			return json_dict

		json_dict = JSONDecoder().decode(string)
		if isinstance(json_dict, list):
			return [Task(**sanity_check(s)) for s in json_dict]
		elif isinstance(json_dict, dict):
			return [Task(**sanity_check(json_dict))]
		else:
			raise IOError('The supplied JSON file does not have a proper structure')

	def write_to_string(self, task_arg):
		## TODO: Make this less hardcoded
		if isinstance(task_arg, (list,tuple,set)):
			d = lambda task: {k:getattr(task,k) for k,dv in Task.__defaults__.items() if getattr(task,k) != dv}
			tasks = [d(t) for t in task_arg]
			return JSONEncoder().encode(tasks)
		elif isinstance(task_arg, Task):
			return JSONEncoder().encode(task_arg)
		else:
			raise TypeError('task_arg is not an iterable containing Task objects or a Task object itself.')


if  __name__ == '__main__':
	jtio = JSONTaskIO()
	tasks = jtio.read_task('resources/generic_models/task_test.json')
	jtio.write_task('resources/generic_models/task_write_test.json', tasks)
