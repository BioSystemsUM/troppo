import numpy as np
from itertools import chain

from cobamp.core.linear_systems import GenericLinearSystem, VAR_CONTINUOUS, VAR_BINARY
from cobamp.core.optimization import LinearSystemOptimizer


class tINIT():
	def __init__(self, S, ub, lb, properties):
		self.S = S
		self.ub = ub
		self.lb = lb
		self.properties = properties
		self.n_metabolites, self.n_reactions = self.S.shape

	def preprocessing(self):
		self.reaction_scores = self.properties['reactions_scores']
		self.present_metabolites = self.properties['present_metabolites']
		self.essential_reactions = self.properties['essential_reactions']
		self.production_weight = self.properties['production_weight']
		self.allow_excretion = self.properties['allow_excretion']
		self.no_reverse_loops = self.properties['no_reverse_loops']

		if self.present_metabolites.size != np.unique(
				self.present_metabolites.size):  # in theory, this should not even happen once, only if by mistake
			print('There are repeated metabolites in the list. We\' apply a unique method and will be used ')
			self.present_metabolites = np.unique(self.present_metabolites)

		self.present_metabolites_unlisted = np.array(
			list(chain(*map(lambda x: [x] if not isinstance(x, list) else x, self.present_metabolites))))

		if self.present_metabolites_unlisted:
			self.metabolites_production = np.ones(self.present_metabolites_unlisted.size) * -1
		# self.metabolites_production[self.present_metabolites] = -1
		else:
			self.metabolites_production = np.array([])
			self.present_metabolites = np.array([])
			self.present_metabolites_unlisted = np.array([])

		# TODO have a function to check if the model is unconstrained
		# TODO change the reversible stuff later to the new to be implemented functions
		self.reversible_reactions = np.where(self.lb < 0)[0]  # change this
		self.essential_reversible_reactions = self.essential_reactions[self.reversible_reactions]
		self.essential_reactions = np.delete(self.essential_reactions, self.essential_reversible_reactions)

		# convert the model to irreversible
		self.irreversible_S, self.irreversible_lb, self.irreversible_ub = self.irreversible_model(
			self.S)  # TODO implement this function
		self.reaction_scores = np.concatenate(self.reaction_scores, self.reaction_scores[self.reversible_reactions])

		if self.no_reverse_loops:
			self.fwd_idx = self.reversible_reactions
			self.bwd_idx = np.concatenate(self.n_reactions,
										  [i + self.n_reactions for i in range(self.irreversible_S.shape[1])])
		else:
			self.fwd_idx = self.essential_reversible_reactions
			self.bwd_idx = [i + self.n_reactions for i in range(len(self.essential_reversible_reactions))]

		# remove the essential reactions from the self.reactions_scores
		self.reaction_scores = np.delete(self.reaction_scores, self.essential_reactions)

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
		self.n_essential_reactions = self.essential_reactions.size
		self.n_non_essential_reactions = self.n_reactions_irrev - self.essential_reactions.size

		# non-essential reactions to produce a fake metabolite
		self.irreversible_S = np.vstack(
			[self.irreversible_S, np.delete(np.eye(self.n_reactions_irrev), self.essential_reactions, 0)])
		# for non-essential, but with a stoichiometry of 1000
		temp = np.eye(self.n_reactions_irrev) * 1000
		new_temp = np.vstack([np.zeros([self.n_metabolites_irrev, self.n_non_essential_reactions]), temp])
		self.irreversible_S = np.hstack(
			[self.irreversible_S, new_temp])

		if self.production_weight != 0:
			temp2 = np.vstack(
				[np.eye(self.n_metabolites_irrev - self.present_metabolites_unlisted.size) * -1,
				 np.zeros([self.n_non_essential_reactions + self.present_metabolites_unlisted.size,
						   self.n_metabolites_irrev - self.present_metabolites_unlisted.size])])

			self.irreversible_S = np.hstack([S, temp2])
			self.n_net_production = self.n_metabolites_irrev - self.present_metabolites_unlisted.size

		else:
			self.n_net_production = 0

		if self.fwd_idx.size > 0:
			self.n_rev_bounds = self.fwd_idx.size
			I = np.eye(self.irreversible_S.shape[1]) * -1
			temp = np.vstack([I[self.fwd_idx,], I[self.bwd_idx,]])
			temp = np.hstack([temp, np.zeros([temp.shape[0], self.irreversible_S.shape[1] - self.n_reactions_irrev])])
			temp = np.hstack([temp, np.eye(self.n_rev_bounds * 2) * 1000])
			temp = np.vstack([temp, np.hstack(
				[np.zeros([self.n_rev_bounds, self.irreversible_S.shape[1]]), np.eye(self.n_rev_bounds) * -1,
				 np.eye(self.n_rev_bounds) * -1])])
			self.irreversible_S = np.hstack(
				[self.irreversible_S, np.zeros([self.irreversible_S.shape[0], self.n_rev_bounds * 2])])
			self.irreversible_S = np.hstack([self.irreversible_S, temp])
		else:
			self.n_rev_bounds = 0

	def build_problem(self):

		# Preparation of the problem to be optimized
		self.problem_blx = np.vstack([self.irreversible_lb, np.zeros(
			self.n_non_essential_reactions + self.n_net_production + (self.n_rev_bounds * 2))])
		# TODO check if this in float; if not, we have to change the dtype of problem_blx, otherwise 0.1 will no exist (will be 0)
		self.problem_blx[self.essential_reactions,] = np.apply_along_axis(lambda x: [i if i > 0.1 else 0.1 for i in x],
																		  1,
																		  self.problem_blx[
																			  self.essential_reactions,])  # this == max(0.1, prob.bl(essentialIndex))

		self.problem_ubx = np.vstack([self.irreversible_ub, np.ones(
			self.n_non_essential_reactions + self.n_net_production + self.n_rev_bounds * 2)])

		self.irreversible_b = np.zeros(self.n_metabolites)

		self.problem_blc = np.vstack(
			[self.irreversible_b, np.ones(self.n_non_essential_reactions), np.zeros(self.n_rev_bounds * 2),
			 np.ones(self.n_rev_bounds) * -1])

		# this part of the code will take care of excretion and reverse loops, if they are true on the parameters
		if self.no_reverse_loops:
			self.rev_ub = np.zeros(self.n_rev_bounds)
			self.rev_ub[self.essential_reversible_reactions] = -1
		else:
			self.rev_ub = np.ones(self.n_rev_bounds) * -1

		if self.allow_excretion:
			self.met_ub = np.array([None] * self.n_metabolites)
		else:
			self.met_ub = np.zeros(self.n_metabolites)

		self.problem_buc = np.vstack(
			[self.met_ub, np.ones(self.n_non_essential_reactions) * 1000, np.ones(self.n_rev_bounds * 2) * 999.9,
			 self.rev_ub])

		self.problem_c = np.vstack([np.zeros(self.n_reactions), self.reaction_scores,
									np.ones(self.n_net_production) * self.production_weight * -1,
									np.zeros(self.n_rev_bounds * 2)])

		self.problem_a = self.irreversible_S

	def solve_problem(self):
		start_index = self.n_metabolites_irrev - self.present_metabolites_unlisted
		problem = GenericLinearSystem(self.problem_a, VAR_CONTINUOUS, self.problem_blx, self.problem_ubx,
									  self.problem_blc, self.problem_buc,
									  ['V' + str(i) for i in range(self.problem_a.shape[1])])
		for i in range(self.present_metabolites_unlisted.size):
			self.problem_blc[start_index + i] = 1
			lso = LinearSystemOptimizer(problem)
			problem.set_constraint_bounds(problem.model.constraints, self.problem_blc, self.problem_buc)
			problem.set_objective(self.problem_c, False)
			solution = lso.optimize()

			if solution.status() != 'optimal':
				self.problem_blc[start_index + i] = 0
			else:
				self.metabolites_production[self.present_metabolites_unlisted[i]] = 1

		allInt = np.hstack([list(range(self.n_reactions + 1, self.n_reactions + self.n_non_essential_reactions)), list(
			range(self.irreversible_S.shape[1] - (self.n_rev_bounds * 2) + 1, self.irreversible_S.shape[1]))])


		problem.set_variable_types([problem.model.variables[i] for i in allInt], VAR_BINARY)
		solution = lso.optimize()

		if solution.status()!='optimal':
			print('The problem is infeasible')
			return # TODO improve




if __name__ == '__main__':
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
