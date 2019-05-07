import pandas as pd
import numpy as np
import re

class Task():
	'''
	This class is to create an object to accommodate the metabolic tasks proposed by this paper
	'''

	def __init__(self, matrix_with_info):
		self._matrix_with_info = matrix_with_info
		self._id = self.parse_id()
		self._system = self.parse_system()
		self._subsystem = self.parse_subsystem()
		self._description = self.parse_description()
		self._fail = self.parse_fail()
		self._input = self.parse_input()
		self._output = self.parse_output()
		self._equation = self.parse_equation()

	@property
	def matrix_with_info(self):
		return self._matrix_with_info

	@matrix_with_info.setter
	def matrix_with_info(self, value):
		self._matrix_with_info = value

	@matrix_with_info.deleter
	def matrix_with_info(self):
		del self._matrix_with_info

	@property
	def id(self):
		return self._id

	@id.setter
	def id(self, value):
		self._id = value

	@id.deleter
	def id(self):
		del self._id

	@property
	def system(self):
		return self._system

	@system.setter
	def system(self, value):
		self._system = value

	@system.deleter
	def system(self):
		del self._system

	@property
	def subsystem(self):
		return self._subsystem

	@subsystem.setter
	def subsystem(self, value):
		self._subsystem = value

	@subsystem.deleter
	def subsystem(self):
		del self._subsystem

	@property
	def description(self):
		return self._description

	@description.setter
	def description(self, value):
		self._description = value

	@description.deleter
	def description(self):
		del self._description

	@property
	def fail(self):
		return self._fail

	@fail.setter
	def fail(self, value):
		self._fail = value

	@fail.deleter
	def fail(self):
		del self._fail

	@property
	def input(self):
		return self._input

	@input.setter
	def input(self, value):
		self._input = value

	@input.deleter
	def input(self):
		del self._input

	@property
	def output(self):
		return self._output

	@output.setter
	def output(self, value):
		self._output = value

	@output.deleter
	def output(self):
		del self._output

	@property
	def equation(self):
		return self._equation

	@equation.setter
	def equation(self, value):
		self._equation = value

	@equation.deleter
	def equation(self):
		del self._equation

	def parse_id(self):
		return int(self._matrix_with_info['ID'].dropna().unique())

	def parse_system(self):
		return self._matrix_with_info['SYSTEM'].dropna().unique()[0]

	def parse_subsystem(self):
		return self._matrix_with_info['SUBSYSTEM'].dropna().unique()[0]

	def parse_description(self):
		return self._matrix_with_info['DESCRIPTION'].dropna().unique()[0]

	def parse_fail(self):
		return int(self._matrix_with_info['SHOULD FAIL'].dropna().unique())

	def parse_input(self):
		return dict(
			zip(self._matrix_with_info['IN'], zip(self._matrix_with_info['IN LB'], self._matrix_with_info['IN UB'])))

	def parse_output(self):
		return dict(
			zip(self._matrix_with_info['OUT'], zip(self._matrix_with_info['OUT LB'], self._matrix_with_info['OUT UB'])))

	def parse_equation(self):
		if self._matrix_with_info['PATHWAYS USED'].dropna().tolist()!=[]:
			base = self._matrix_with_info['PATHWAYS USED'].dropna().tolist()[0]
			first_split = re.split( " \W*- *(?=[A-Z])", base)
			dict_equations= {s[0:s.find(' (')]:s[s.find(' (')+2:-1] for s in first_split}
			return dict_equations
		return {}

if __name__ == '__main__':
	path = './data/Nathan2019ConsensusPaper/'
	z = pd.read_excel(path + 'pcbi.1006867.s005.xlsx')
	idx = np.where(z.ID.isna() == False)[0]

	t = Task(z.iloc[idx[0]:idx[1],:])