import abc
from io import TextIOBase, StringIO
from json import JSONDecoder, JSONEncoder
import warnings
from troppo.tasks.core import Task
import pandas as pd

class TaskIO(object):

	def read_task(self, buffer_or_path, binary_mode=False):
		if isinstance(buffer_or_path, TextIOBase):
			return self.read_from_string(buffer_or_path.read())

		elif isinstance(buffer_or_path, str):
			with open(buffer_or_path, 'r'+('b' if binary_mode else '')) as f:
				return self.read_from_string(f.read())

		else:
			raise TypeError('Invalid buffer or path')

	def write_task(self, buffer_or_path, task_arg, binary_mode=False):
		task_string = self.write_to_string(task_arg)
		if isinstance(buffer_or_path, TextIOBase):
			buffer_or_path.write(task_string)
		elif isinstance(buffer_or_path, str):
			with open(buffer_or_path, 'w'('b' if binary_mode else '')) as f:
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

class ExcelTaskIO(TaskIO):
	def read_task(self, buffer_or_path):
		return super().read_task(buffer_or_path, True)

	def read_from_string(self, string):
		shd_fail = ['SHOULD FAIL']
		bounds = ['LB', 'UB']
		inflows = [('IN ' + x).strip() for x in [''] + bounds]
		outflows = [('OUT ' + x).strip() for x in [''] + bounds]
		eqs = [('EQU ' + x).strip() for x in [''] + bounds]
		flx_const = [('CHANGED ' + x).strip() for x in ['RXN'] + bounds]

		core_info = shd_fail + inflows + outflows + eqs + flx_const
		tdf = pd.read_excel((string))
		valid_tdf = tdf.loc[tdf.iloc[:, 0] != '#', :].iloc[:, 1:]

		# create unique ids
		from collections import Counter
		id_counter = Counter()
		lbs, ubs = [[inflows[k], outflows[k], eqs[k]] for k in [1, 2]]
		valid_tdf[lbs] = valid_tdf[lbs].fillna(0)
		valid_tdf[ubs] = valid_tdf[ubs].fillna(1000)

		real_id = []
		for i in valid_tdf['ID'].fillna(-1):
			if i != -1:
				id_counter[i] += 1
				real_id.append(str(i) + str(id_counter[i]))
			else:
				real_id.append(real_id[-1])
		valid_tdf['ID'] = real_id

		task_info = {}
		for task_id, rows in tuple(valid_tdf.groupby('ID')):
			annotations = rows.drop(columns=['ID'] + core_info).fillna('') \
				.apply(lambda x: ','.join([k for k in x if k != ''])).to_dict()
			task_info[task_id] = {'annotations': annotations}
			task_info[task_id]['name'] = task_id
			for dtname, datatype in zip(['inflow_dict', 'outflow_dict', 'reaction_dict', 'flux_constraints'],
			                            [inflows, outflows, eqs, flx_const]):
				task_info[task_id][dtname] = {}
				for i, row in rows[datatype].iterrows():
					if not row.isnull().iloc[0]:
						task_info[task_id][dtname].update({k: [row[1], row[2]] for k in row.iloc[0].split(';')})
			task_info[task_id] = {k: v for k, v in task_info[task_id].items() if len(v) > 0}
			# parse reactions
			if 'reaction_dict' in task_info[task_id].keys():
				new_reactions_dict = {}
				for rxstr, bounds in task_info[task_id]['reaction_dict'].items():
					if '<=>' in rxstr:
						delim, rev = '<=>', True
					else:
						delim, rev = '=>', False
					rctns, prods = map(lambda x: x.strip(), rxstr.split(delim))
					coef_dict = {}
					for l, z in zip([rctns, prods], [-1, 1]):
						coef_dict.update(
							{k[0]: (1 if len(k) == 1 else float(k[1])) * z for k in
							 [y.strip().split(' ') for y in l.split(' + ')]})
					new_reactions_dict['reaction_' + task_id + '_' + str(len(new_reactions_dict) + 1)] = (
					coef_dict, bounds)
				task_info[task_id]['reaction_dict'] = new_reactions_dict
			task_info[task_id]['should_fail'] = rows[shd_fail[0]].fillna(False).iloc[0]

		return [Task(**v) for k, v in task_info.items()]

	def write_to_string(self, task):
		raise Exception('Not yet implemented!')

if  __name__ == '__main__':
	# jtio = JSONTaskIO()
	# tasks = jtio.read_task('resources/generic_models/task_test.json')
	# jtio.write_task('resources/generic_models/task_write_test.json', tasks)

	etio = ExcelTaskIO()
	tasks = etio.read_task('shared/task_sets/metabolicTasks_Essential.xlsx')
	import pandas as pd