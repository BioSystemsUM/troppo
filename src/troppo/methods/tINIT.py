import numpy as np
import scipy.sparse as sprs
from numpy import int_
from itertools import chain

from cobamp.core.linear_systems import GenericLinearSystem, VAR_CONTINUOUS, VAR_BINARY
from cobamp.core.optimization import LinearSystemOptimizer

from cobamp.core.models import make_irreversible_model, make_irreversible_model_raven
from troppo.methods.base import ContextSpecificModelReconstructionAlgorithm, PropertiesReconstruction
from troppo.reconstruction_properties import is_list_else_empty, if_none_return_list

class tINITProperties(PropertiesReconstruction):
	def __init__(self, reactions_scores, present_metabolites=[], essential_reactions=[], production_weight=0.5,
				 allow_excretion=False,
				 no_reverse_loops=False, solver=None):
		# TODO check later if the params input might be necessary to include for the troppo
		new_mandatory = {
			'solver': lambda x: isinstance(x, str)
		}
		new_optional = {
			'reactions_scores': lambda x: is_list_else_empty(x),
			'present_metabolites': lambda x: is_list_else_empty(x),
			'essential_reactions': lambda x: is_list_else_empty(x),
			'production_weight': lambda x: isinstance(x, float),
			'allow_excretion': lambda x: isinstance(x, bool),
			'no_reverse_loops': lambda x: isinstance(x, bool),
		}

		super().__init__()

		self.add_new_properties(new_mandatory, new_optional)

		properties_dict_list = ["reactions_scores", "present_metabolites", "essential_reactions", "production_weight",
								"allow_excretion", "no_reverse_loops", "solver"]
		properties_list = [reactions_scores, if_none_return_list(present_metabolites),
						   if_none_return_list(essential_reactions), if_none_return_list(production_weight),
						   allow_excretion, no_reverse_loops, solver]
		prop_dict = dict(zip(properties_dict_list, properties_list))

		[self.add_if_not_none(*k) for k in prop_dict.items()]

	@staticmethod
	def from_integrated_scores(scores, **kwargs):
		return tINITProperties(reactions_scores=scores, **kwargs)


class tINIT(ContextSpecificModelReconstructionAlgorithm):
	properties_class = tINITProperties
	def __init__(self, S, lb, ub, properties):
		super().__init__(S, lb, ub, properties)
		self.S = np.array(S)
		self.lb = np.array(lb)
		self.ub = np.array(ub)
		self.properties = properties
		self.n_metabolites, self.n_reactions = self.S.shape

	def preprocessing(self):
		self.reaction_scores = np.array(self.properties['reactions_scores'])
		self.present_metabolites = np.array(self.properties['present_metabolites'])

		self.essential_reactions_idx = np.array(self.properties['essential_reactions'])
		self.essential_reactions = np.zeros(self.n_reactions).astype(bool)
		if self.essential_reactions_idx.size != 0:
			self.essential_reactions[self.essential_reactions_idx] = True

		self.production_weight = self.properties['production_weight']
		self.allow_excretion = self.properties['allow_excretion']
		self.no_reverse_loops = self.properties['no_reverse_loops']

		if self.present_metabolites.size != np.unique(
				self.present_metabolites.size):  # in theory, this should not even happen once, only if by mistake
			print('There are repeated metabolites in the list. We\' apply a unique method and will be used ')
			self.present_metabolites = np.unique(self.present_metabolites)

		self.present_metabolites_unlisted = np.array(
			list(chain(*map(lambda x: [x] if not isinstance(x, list) else x, self.present_metabolites))))

		if self.present_metabolites_unlisted.size > 0:
			self.metabolites_production = np.ones(self.present_metabolites_unlisted.size) * -1
		# self.metabolites_production[self.present_metabolites] = -1
		else:
			self.metabolites_production = np.array([])
			self.present_metabolites = np.array([])
			self.present_metabolites_unlisted = np.array([])

		# TODO have a function to check if the model is unconstrained
		# TODO change the reversible stuff later to the new to be implemented functions
		self.reversible_reactions = np.where(self.lb < 0)[0] # change this TODO should this also have self.ub > 0?
		# self.reversible_reactions = np.where((np.logical_and(self.lb < 0, self.ub>0)))[0]
		self.essential_reversible_reactions = self.essential_reactions_idx[
			np.isin(self.essential_reactions_idx, self.reversible_reactions)]
		self.essential_reactions_idx = np.setdiff1d(self.essential_reactions_idx, self.essential_reversible_reactions)

		# convert the model to irreversible
		self.irreversible_S, self.irreversible_lb, self.irreversible_ub, self.rev_map = make_irreversible_model_raven(self.S,
																												self.lb,
																												self.ub, False)

		self.reaction_scores = np.hstack([self.reaction_scores, self.reaction_scores[
			self.reversible_reactions]])

		if self.no_reverse_loops:
			self.fwd_idx = self.reversible_reactions
			self.bwd_idx = np.array([self.rev_map[i][1] for i in self.reversible_reactions]) #TODO check this part of no reverse loops

		else:
			self.fwd_idx = self.essential_reversible_reactions
			# self.bwd_idx = np.array([i + (self.n_reactions) for i in self.essential_reversible_reactions])
			self.bwd_idx = np.array([self.rev_map[i][1] for i in self.essential_reversible_reactions])

		# remove the essential reactions from the self.reactions_scores
		self.reaction_scores = np.delete(self.reaction_scores, self.essential_reactions_idx)

		# this part will create fake metabolites for the self.present_metabolites to insert on the self.irreversible_S
		# TODO this should also be a function in a superclass
		if self.present_metabolites_unlisted.size > 0:
			metabolites_to_add = [i + self.irreversible_S.shape[0] for i in
								  range(self.present_metabolites_unlisted.size)]

		pre_shape_irreversible_iterable = self.irreversible_S.shape[0]

		self.irreversible_S = np.vstack(
			[self.irreversible_S, np.zeros([self.present_metabolites_unlisted.size, self.irreversible_S.shape[1]])])

		# change the values for the new added metabolites
		# TODO check this for to make sure that the unlisted present metabolites have the same location as the listed present metabolites
		for metabolites in self.present_metabolites:
			if type(metabolites) in [list, tuple]:
				# check the reactions for all the metabolites in the list/tuple
				temp_S = self.S[metabolites,]
				temp_idx, temp_stoi = np.nonzero(temp_S), temp_S[np.nonzero(temp_S)]
				unique_reactions = np.unique(temp_idx[1])
				unique_stoi = [np.sum(temp_S[:, i]) for i in range(unique_reactions.size)]

				for metabolite in metabolites:
					self.irreversible_S[pre_shape_irreversible_iterable, temp_idx] = unique_stoi
					pre_shape_irreversible_iterable = pre_shape_irreversible_iterable + 1
			# temp_values = self.S[metabolite,]
			# temp_idx, temp_stoi = np.nonzero(temp_values), temp_values[np.nonzero(temp_values)]
			else:
				temp_values = self.S[metabolites,]
				temp_idx, temp_stoi = np.nonzero(temp_values), temp_values[np.nonzero(temp_values)]
				self.irreversible_S[pre_shape_irreversible_iterable, temp_idx] = temp_stoi
				pre_shape_irreversible_iterable = pre_shape_irreversible_iterable + 1

		# TODO CHECK THIS LATER TO SEE IF THIS DOESN'T FUCKING MESS UP

		# TODO for this block, check if the values are updated throughout the modifications, but I do not think so and its the way it should be
		# useful numbers
		self.n_metabolites_irrev, self.n_reactions_irrev = self.irreversible_S.shape
		self.n_essential_reactions = self.essential_reactions_idx.size
		self.n_non_essential_reactions = self.n_reactions_irrev - self.n_essential_reactions
		self.non_essential = np.setdiff1d(list(range(self.irreversible_S.shape[1])), self.essential_reactions_idx)

		self.v_vector_reactions = list(self.rev_map.keys())
		for rev_r in self.reversible_reactions: self.v_vector_reactions.append(str(rev_r) + '_r')

		revs = np.setdiff1d(self.non_essential, list(self.rev_map.keys()))

		# revs_names = []
		# for r in range(revs.size):
		# 	for k, v in self.rev_map.items():
		# 		if isinstance(v, (tuple, list)):
		# 			if revs[r] in v:
		# 				revs_names.append(str(k) + '_r')

		revs_names = [str(k)+'_r' for k,v in self.rev_map.items() if isinstance(v, (tuple, list)) and v[1] in revs]

		self.non_essential = self.non_essential.astype(str)
		revs = revs.astype(str)

		for r in range(revs.size):
			temp = np.where(self.non_essential == revs[r])[0]
			self.non_essential[temp] = revs_names[r]

		## TODO: The blocks below this will need a rework as larger models will use too much memory
		# non-essential reactions to produce a fake metabolite

		ids_to_keep = np.setdiff1d(np.arange(self.n_reactions_irrev), self.essential_reactions_idx)
		matrix_to_add = sprs.eye(self.n_reactions_irrev, format='csc')[:, ids_to_keep]

		self.irreversible_S = sprs.vstack(
			[self.irreversible_S, matrix_to_add])

		# for non-essential, but with a stoichiometry of 1000
		temp = sprs.eye(self.n_non_essential_reactions, format='csc') * 1000
		new_temp = sprs.vstack([np.zeros([self.n_metabolites_irrev, self.n_non_essential_reactions]), temp])
		self.irreversible_S = sprs.hstack(
			[self.irreversible_S, new_temp])

		if self.production_weight != 0:
			zmat = sprs.csc_matrix(np.zeros([self.n_non_essential_reactions + self.present_metabolites_unlisted.size,
					  self.n_metabolites_irrev - self.present_metabolites_unlisted.size]))
			temp2 = sprs.vstack(
				[sprs.eye(self.n_metabolites_irrev - self.present_metabolites_unlisted.size, format='csc') * -1,
				 zmat])

			self.irreversible_S = sprs.hstack([self.irreversible_S, temp2])
			self.n_net_production = self.n_metabolites_irrev - self.present_metabolites_unlisted.size
			self.net_production = np.setdiff1d(list(range(self.n_metabolites_irrev)), self.present_metabolites_unlisted)

		else:
			self.net_production = []
			self.n_net_production = 0

		if self.fwd_idx.size > 0:
			self.n_rev_bounds = self.fwd_idx.size
			I = sprs.eye(self.n_reactions_irrev) * -1
			temp = sprs.vstack([I[self.fwd_idx,], I[self.bwd_idx,]])
			# padding
			temp = sprs.hstack([temp, sprs.csc_matrix(np.zeros([temp.shape[0], self.irreversible_S.shape[1] - self.n_reactions_irrev]))])
			temp = sprs.hstack([temp, sprs.eye(self.n_rev_bounds * 2, format='csc') * 1000])
			temp = sprs.vstack([temp, sprs.hstack(
				[sprs.csc_matrix(np.zeros([self.n_rev_bounds, self.irreversible_S.shape[1]])), sprs.eye(self.n_rev_bounds, format='csc') * -1,
				 sprs.eye(self.n_rev_bounds, format='csc') * -1])])
			self.irreversible_S = sprs.hstack(
				[self.irreversible_S, sprs.csc_matrix(np.zeros([self.irreversible_S.shape[0], self.n_rev_bounds * 2]))])
			self.irreversible_S = sprs.vstack([self.irreversible_S, temp])
		else:
			self.n_rev_bounds = 0

	def build_problem(self):

		# Preparation of the problem to be optimized
		self.problem_blx = sprs.hstack([self.irreversible_lb, sprs.csc_matrix(np.zeros(
			self.n_non_essential_reactions + self.n_net_production + (self.n_rev_bounds * 2)))])
		# TODO check if this in float; if not, we have to change the dtype of problem_blx, otherwise 0.1 will no exist (will be 0)
		if self.essential_reactions_idx.size > 0:
			self.problem_blx[self.essential_reactions_idx] = np.array(
				list(map(lambda x: x if x > 0.1 else 0.1, self.problem_blx[self.essential_reactions_idx])))

		self.problem_bux = sprs.hstack([self.irreversible_ub, np.ones(
			self.n_non_essential_reactions + self.n_net_production + self.n_rev_bounds * 2)])

		self.irreversible_b = np.zeros(self.n_metabolites_irrev)

		self.problem_blc = np.hstack(
			[self.irreversible_b, np.ones(self.n_non_essential_reactions), np.zeros(self.n_rev_bounds * 2),
			 np.ones(self.n_rev_bounds) * -1])

		# this part of the code will take care of excretion and reverse loops, if they are true on the parameters
		if self.no_reverse_loops:
			self.rev_ub = np.zeros(self.n_rev_bounds)
			# TODO this is just a fix; REALLY FIX THIS LATER
			for i in range(self.essential_reversible_reactions.size):
				self.rev_ub[i] = -1
		else:
			self.rev_ub = np.ones(self.n_rev_bounds) * -1

		if self.allow_excretion:
			self.met_ub = np.array([None] * self.n_metabolites_irrev)
		else:
			self.met_ub = np.zeros(self.n_metabolites_irrev)

		self.problem_buc = np.hstack(
			[self.met_ub, np.ones(self.n_non_essential_reactions) * 1000, np.ones(self.n_rev_bounds * 2) * 999.9,
			 self.rev_ub])

		self.problem_c = np.hstack([np.zeros(self.n_reactions_irrev), self.reaction_scores,
									np.ones(self.n_net_production) * self.production_weight * -1,
									np.zeros(self.n_rev_bounds * 2)])

		self.problem_a = self.irreversible_S

	def solve_problem(self):

		# len_reactions = [len(self.irreversible_lb), self.n_non_essential_reactions, self.n_net_production,
		# 				 self.n_rev_bounds, self.n_rev_bounds]
		# prefix = ['v_', 'ne_', 'net_prod_', 'revF_', 'revB_']
		#
		# rx_names_problem2 = list(chain(*[[k + str(i) for i in range(v)] for k, v in zip(prefix, len_reactions)]))

		len_reactions = [np.array(self.v_vector_reactions), self.non_essential, self.net_production,
						 self.fwd_idx, self.fwd_idx]
		prefix = ['v_', 'ne_', 'net_prod_', 'revF_', 'revB_']

		rx_names_problem = list(chain(*[[k + str(i) for i in v] for k, v in zip(prefix, len_reactions)]))

		start_index = self.n_metabolites_irrev - self.present_metabolites_unlisted.size
		problem = GenericLinearSystem(self.problem_a.toarray(), VAR_CONTINUOUS, self.problem_blx, self.problem_bux,
									  self.problem_blc, self.problem_buc,
									  rx_names_problem, self.properties['solver'])

		## TODO: Remove this later
		if self.properties['solver'] == 'CPLEX':
			problem.model.problem.parameters.mip.tolerances.mipgap.set(1e-9)
		elif self.properties['solver'] == 'GUROBI':
			problem.model.problem.Params.MIPGap = 1e-9

		problem.model.configuration.tolerances.feasibility = 1e-8
		problem.model.configuration.tolerances.optimality = 1e-8
		problem.model.configuration.verbosity = 3

		lso = LinearSystemOptimizer(problem)
		problem.write_to_lp('tINIT_test.lp')
		if self.present_metabolites_unlisted.size > 0:
			for i in range(self.present_metabolites_unlisted.size):
				self.problem_buc[start_index + i] = 1
				problem.set_constraint_bounds(problem.model.constraints, self.problem_blc, self.problem_buc)
				problem.set_objective(self.problem_c, True)
				solution = lso.optimize()

				if solution.status() != 'optimal':
					self.problem_buc[start_index + i] = 0
				else:
					self.metabolites_production[i] = 1

		allInt = np.hstack([list(range(self.n_reactions_irrev, self.n_reactions_irrev + self.n_non_essential_reactions)), list(
			range(self.irreversible_S.shape[1] - (self.n_rev_bounds * 2) + 1, self.irreversible_S.shape[1]))])

		# return lso.optimize()

		problem.set_variable_types([problem.model.variables[int_(i)] for i in allInt], VAR_BINARY)
		problem.set_objective(self.problem_c, True)
		problem.write_to_lp('tINIT_test.lp')
		# print(problem.model.to_lp())
		solution = lso.optimize()
		print(solution.objective_value())

		if solution.status() != 'optimal':
			print('The problem is infeasible')
			# return solution  # TODO improve
			return
		# return solution

		# used_reactions = np.array([])

		# for k, v in solution.var_values().items(): print(k, v)

		used_reactions = np.array([int_(k.split('_')[1]) for k in list(solution.var_values().keys()) if 'ne_' in k and solution.var_values()[k] < 0.1])

		# for k, v in solution.var_values().items():
		# 	if 'ne_' in k:
		# 		# if '_r' in k:
		# 		# if v < 1 - 1e-6:
		# 		if v < 0.1:
		# 			used_reactions = np.append(used_reactions, int_(k.split('_')[1]))

		return np.append(used_reactions, self.essential_reactions_idx)

	def run_tINIT(self):
		self.preprocessing()
		self.build_problem()
		res = self.solve_problem()
		return np.unique(np.int_(np.sort(res)))

	def run(self):
		return self.run_tINIT()



if __name__ == '__main__':
	import numpy as np

	S = np.array([[1, -1, 0, 0, -1, 0, -1, 0, 0],
				  [0, 1, -1, 0, 0, 0, 0, 0, 0],
				  [0, 1, 0, 1, -1, 0, 0, 0, 0],
				  [0, 0, 0, 0, 0, 1, -1, 0, 0],
				  [0, 0, 0, 0, 0, 0, 1, -1, 0],
				  [0, 0, 0, 0, 1, 0, 0, 1, -1]])

	a = S[2,]
	n_a_idx, n_a = np.nonzero(a), a[np.nonzero(a)]
	x = S[[1, 2, 3],]
	n_x_idx, n_x = np.nonzero(x), x[np.nonzero(x)]
	unique_reactions = np.unique(n_x_idx[1])
	z = np.array([np.sum(x[:, i]) for i in range(unique_reactions.size)], dtype=float)
	z = np.vstack([z, z, z, z, z, z])
	z[[0, 2, 4],] = np.apply_along_axis(lambda x: [i if i > 0.1 else 0.1 for i in x], 1, z[[0, 2, 4],])
	l = (lambda x: [i if i > 0.1 else 0.1 for i in x])
	l([0, 1, 2])


