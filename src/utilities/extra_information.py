import pandas as pd
import numpy as np
from utilities.task import Task
import requests, zipfile, io


class InformationFromSeveralArticles:

    def __init__(self):
        pass

    @staticmethod
    def hmr_to_recon_reaction_converter():
        """
        This function maps the HMR reactions to the Recon reactions.
        Returns
        -------
        mapping: pandas.DataFrame
            Dataframe with the mapping.
        """
        hmr_zip = requests.get("http://www.metabolicatlas.org/assets/hmr/HMRdatabase2_00.xlsx"
                               "-07f0f2b5f71a03b0d6efdea752a9cebb.zip")
        with zipfile.ZipFile(io.BytesIO(hmr_zip.content)) as my_zip:
            hmr = pd.ExcelFile(my_zip.open('HMRdatabase2_00.xlsx'))

        rx_hmr = hmr.parse(['RXNS'])['RXNS']
        mapping = rx_hmr[['RXNID', 'BIGG DATABASE ID', 'COMPARTMENT']]
        return mapping

    def load_medium_rpmi1640(self, path_rpmi='./data/extra_information/rpmi1640.tsv'):
        """
        This function loads the rpmi-1640 information from the paper https://www.nature.com/articles/srep45557
        NOTE: this information is present in the Supplementary Information, but the data is present on a .doc, so we
        had to get that information, put it on a .tsv and then save it.
        Parameters:
        -----------
        path_rpmi: str
            Path to the rpmi-1640 file.
        Returns:
        --------
        rpmi: pandas.DataFrame
            Dataframe with the rpmi-1640 information.
        """

        rpmi = pd.read_csv(path_rpmi, delimiter='\t').set_index('Reaction ID', drop=True)
        reaction_map = self.hmr_to_recon_reaction_converter().set_index('RXNID', drop=True)
        rpmi[reaction_map.columns] = reaction_map.loc[rpmi.index.values, :]

        return rpmi


class Nathan2019ConsensusPaper:
    """
    This class has the purpose to import and make accessible the information present in the paper and in the
    supplementary material.
    LINK : https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006867
    """

    def __init__(self, root_path='projects/breast_mcf7/data/nl_consensus/'):
        self.root_path = root_path

    def load_medium_constraints(self):
        """
        This function loads the medium constraints from the paper.

        Returns
        -------
        final: pandas.DataFrame
            Dataframe with the medium constraints.
        """
        original = pd.read_excel(self.root_path + 'pcbi.1006867.s003.xlsx', index_col=0)
        cell_lines = [col for col in original.columns if 'Unnamed' not in col]
        gen_names = lambda x: (str(x) + '_LB', str(x) + '_UB')
        new_cols = [bound for pair in [gen_names(cell) for cell in cell_lines] for bound in pair]
        original.columns = new_cols
        final = original.drop(np.nan)
        return final

    def load_biomass(self):
        """
        Function to load the biomass equations for iHsa and Recon2.2 models (based on experimental data)
        and cell lines growth rate.

        Returns
        -------
        ihsa: pandas.DataFrame
            Biomass equations for iHsa.
        recon2_2: pandas.DataFrame
            Biomass equations for Recon2.2.
        real_biomass: pandas.DataFrame
            Cell lines growth rates.
        """
        original = pd.read_excel(self.root_path + 'pcbi.1006867.s004.xlsx', header=1)
        ihsa, recon2_2, real_biomass = original.iloc[:, [0, 1]].dropna(), original.iloc[:, [3, 4]], \
            original.iloc[:, [6, 7]].dropna()
        return ihsa, recon2_2, real_biomass

    def load_metabolic_tasks(self):
        """
        This function loads the metabolic tasks from the paper.

        Returns
        -------
        dict
            Dictionary with the metabolic tasks.
        """
        original = pd.read_excel(self.root_path + 'pcbi.1006867.s005.xlsx')
        idx = np.where(original.ID.isna() == False)[0]
        idx = np.append(idx, original.shape[0] - 1)
        return {i + 1: Task(original.iloc[idx[i]:idx[i + 1], :]) for i in range(len(idx) - 1)}

    def load_hallmarks_genes(self):
        """
        This function loads the hallmark genes from the paper.

        Returns
        -------
        original: pandas.DataFrame
            Dataframe with the hallmark genes.
        """
        original = pd.read_excel(self.root_path + 'pcbi.1006867.s013.xlsx', header=1)
        return original


if __name__ == '__main__':
    data_path = './data/Nathan2019ConsensusPaper/'

    n = Nathan2019ConsensusPaper(data_path)
    x = n.load_medium_constraints()
    c, v, b = n.load_biomass()
