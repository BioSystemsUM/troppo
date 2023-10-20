from typing import Iterable

from cobamp.core.linear_systems import get_default_solver
from cobamp.algorithms.kshortest import *
from troppo.methods.base import GapfillAlgorithm, GapfillProperties
from numpy import ndarray

DEFAULT_CONFIG = KShortestProperties()
DEFAULT_CONFIG[K_SHORTEST_MPROPERTY_METHOD] = K_SHORTEST_METHOD_ITERATE
DEFAULT_CONFIG[K_SHORTEST_OPROPERTY_BIG_M_CONSTRAINTS] = True
DEFAULT_CONFIG[K_SHORTEST_OPROPERTY_BIG_M_VALUE] = 1e6
DEFAULT_CONFIG[K_SHORTEST_OPROPERTY_FORCE_NON_CANCELLATION] = True
DEFAULT_CONFIG[K_SHORTEST_OPROPERTY_MAXSOLUTIONS] = 1


class EFMGapfillProperties(GapfillProperties):
    """
    Properties for the EFM Gap-filling algorithm.

    Parameters:
    -----------
    avbl_fluxes: list
        List of available fluxes.
    lsystem_args: dict
        Dictionary of arguments to be passed to the IrreversibleLinearPatternSystem.
    solver: str
        Solver to be used by the IrreversibleLinearPatternSystem.
    kshproperties: KShortestProperties
        Properties for the KShortestEFMAlgorithm.
    """

    def __init__(self, avbl_fluxes: list, lsystem_args: dict, solver: str = get_default_solver(),
                 kshproperties: KShortestProperties = DEFAULT_CONFIG):
        super().__init__()
        self.add_new_properties({'kshproperties': KShortestProperties}, {})
        self['avbl_fluxes'] = avbl_fluxes
        self['lsystem_args'] = lsystem_args
        self['solver'] = solver
        self['kshproperties'] = kshproperties


class EFMGapfill(GapfillAlgorithm):
    """
    Gap-filling algorithm based on the KShortestEFMAlgorithm.

    Parameters:
    -----------
    S: ndarray
        Stoichiometric matrix.
    lb: ndarray
        Lower bounds.
    ub: ndarray
        Upper bounds.
    properties: EFMGapfillProperties
        Properties for the algorithm.
    """
    properties_class = EFMGapfillProperties

    def __init__(self, S: ndarray, lb: ndarray, ub: ndarray, properties: EFMGapfillProperties):
        super().__init__(S, lb, ub, properties)
        self.__S, self.__lb, self.__ub = S, lb, ub
        self.properties = properties
        self.properties['kshproperties'][K_SHORTEST_MPROPERTY_TYPE_EFP] = False

    def get_enumerator(self, S: ndarray, lb: ndarray, ub: ndarray, avbl_fluxes: list, solver: str,
                       lsystem_args: dict) -> Iterable:
        """
        Get the enumerator for the EFM Gap-filling algorithm.

        Parameters
        ----------
        S: ndarray
            Stoichiometric matrix.
        lb: ndarray
            Lower bounds.
        ub: ndarray
            Upper bounds.
        avbl_fluxes: list
            Available fluxes.
        solver: str
            Solver to be used by the IrreversibleLinearPatternSystem.
        lsystem_args: dict
            Dictionary of arguments to be passed to the IrreversibleLinearPatternSystem.

        Returns
        -------
        enumerator: Iterable
            Enumerator for the EFM Gap-filling algorithm.
        """
        ils = IrreversibleLinearPatternSystem(S, lb, ub, subset=avbl_fluxes, solver=solver, **lsystem_args)
        ksefm = KShortestEFMAlgorithm(self.properties['kshproperties'])
        enumerator = ksefm.get_enumerator(ils, [], [], True)
        return enumerator

    def gapfill(self, avbl_fluxes: list, lsystem_args: dict, solver: str) -> Iterable:
        """
        Gap-filling algorithm based on the KShortestEFMAlgorithm.

        Parameters
        ----------
        avbl_fluxes: list
            Available fluxes.
        lsystem_args: dict
            Dictionary of arguments to be passed to the IrreversibleLinearPatternSystem.
        solver: str
            Solver to be used by the IrreversibleLinearPatternSystem.

        Returns
        -------
        enumerator: Iterable
            Enumerator for the EFM Gap-filling algorithm.
        """
        return self.get_enumerator(self.__S, self.__lb, self.__ub, avbl_fluxes, solver, lsystem_args)

    def run(self) -> list:
        """
        Run the EFM Gap-filling algorithm.

        Returns
        -------
        result: list
            Indices of active indicator variables (maps with variables on the original stoichiometric matrix)
        """
        enm = self.gapfill(self.properties['avbl_fluxes'],
                           lsystem_args=self.properties['lsystem_args'],
                           solver=self.properties['solver'])
        if self.properties['kshproperties'][K_SHORTEST_MPROPERTY_METHOD] == K_SHORTEST_METHOD_ITERATE:
            slist = list(enm)
        elif self.properties['kshproperties'][K_SHORTEST_MPROPERTY_METHOD] == K_SHORTEST_METHOD_POPULATE:
            slist = list(chain(*enm))
        else:
            warnings.warn('Could not find any solution')
            slist = []
        result = [s.get_active_indicator_varids() for s in slist]

        return result


if __name__ == '__main__':
    from urllib.request import urlretrieve
    from cobra.io import read_sbml_model
    from cobra.core import Reaction, Metabolite
    from cobamp.wrappers.cobra import COBRAModelObjectReader
    from random import sample, randint
    from troppo.methods_wrappers import GapfillWrapper

    BIOMASS_RX_NAME = 'BIOMASS_Ecoli_core_w_GAM'
    BIOMASS_TRANS_NAME = 'BIOMASSt'
    BIOMASS_DRAIN_NAME = 'EX_biomass_e'
    BIOMASS_CYT_NAME = 'biomass_c'
    BIOMASS_EXT_NAME = 'biomass_e'


    def random_knockouts(model, biomass_reactions, drains):
        exclude = biomass_reactions | drains
        choices = set(r.id for r in model.reactions) - set(exclude)
        missing_set = set(sample(choices, k=randint(len(choices) // 4, (3 * len(choices)) // 4)))
        return missing_set


    def simulate_model_with_kos(model, knockouts):
        with model as temp_model:
            for rx in knockouts:
                temp_model.reactions.get_by_id(rx).knock_out()
            solution = temp_model.optimize()
        return solution


    def get_ecoli_core_model():
        model_url = 'http://bigg.ucsd.edu/static/models/e_coli_core.xml'
        model_path, _ = urlretrieve(model_url)
        model = read_sbml_model(model_path)
        b_c, b_e = (Metabolite(BIOMASS_CYT_NAME, name='biomass (cytosol)'),
                    Metabolite(BIOMASS_EXT_NAME, name='biomass (extracellular)'))

        model.reactions.get_by_id(BIOMASS_RX_NAME).add_metabolites({b_c: 1})

        b_trans = Reaction(BIOMASS_TRANS_NAME, name='Biomass transport')
        b_trans.add_metabolites({b_c: -1, b_e: 1})

        b_drain = Reaction(BIOMASS_DRAIN_NAME, name='Biomass drain')
        b_drain.add_metabolites({b_e: -1})

        model.add_reactions([b_trans, b_drain])

        return model


    ec_model = get_ecoli_core_model()

    objreader = COBRAModelObjectReader(ec_model)
    biomass_rx_ids = {BIOMASS_RX_NAME, BIOMASS_TRANS_NAME, BIOMASS_DRAIN_NAME}
    non_consumed = {b.id for b in ec_model.boundary if b.bounds[0] >= 0}
    drain_rx_ids = {r.id for r in ec_model.boundary}
    consumed = {BIOMASS_EXT_NAME}

    cobamp_model = objreader.to_cobamp_cbm('CPLEX')
    gpfl_wrapper = GapfillWrapper(model=ec_model)
    failed = []

    for i in range(100):
        kos = random_knockouts(ec_model, biomass_rx_ids, drain_rx_ids)
        ls_override = {'produced': consumed,
                       'non_consumed': [list(ec_model.reactions.get_by_id(r).metabolites)[0].id for r in non_consumed]}
        sol = set(gpfl_wrapper.run(avbl_fluxes=kos, ls_override=ls_override, algorithm='efm')[0])

        with cobamp_model as m:
            for r in (kos - sol):
                m.set_reaction_bounds(r, lb=0, ub=0)
            sol_after = m.optimize({BIOMASS_DRAIN_NAME: 1})
            for r in sol:
                m.set_reaction_bounds(r, lb=0, ub=0)
            sol_before = m.optimize({BIOMASS_DRAIN_NAME: 1})

        if (not sol_after.objective_value() > 0) or sol_after.status() == 'infeasible':
            failed.append((kos, sol, sol_before, sol_after, sol))
