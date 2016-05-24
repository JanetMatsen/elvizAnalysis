import sys
print(sys.version)

module_dir = "../"
sys.path.append(module_dir)

import unittest

import pandas as pd

import abundance_plot_utils

def load_animal_data():
    return pd.read_csv("./summarised_animals.txt", sep='\t')

class testDeleteRowsForTaxa(unittest.TestCase):
    def test_Orca_Deletion(self):
        taxa_level = "Species"
        taxa_name = "Orcinus orca"
        # Load the data
        animal_df = pd.read_csv("./summarised_animals.txt", sep='\t')
        self.assertTrue(taxa_name in animal_df[taxa_level].unique())

        # Now remove that taxa and make sure it is gone.
        animal_df = abundance_plot_utils.delete_rows_for_taxa(
        animals, 'Species', 'Orcinus orca')
        self.assertFalse(taxa_name in animal_df[taxa_level].unique())
        print("unique values left: {}".format(animal_df[taxa_level].unique()))


class testAggregateMixedTaxonomy(unittest.TestCase):
    def testLumpingIntoOther(self):
        # TODO: incomplete. But there is a check in aggregate_mixed_taxonomy():
        # `assert (sample_sums < 1.001).all()`
        pass


if __name__ == "__main__":
    # add the main dir to the path
    print('run unit test for ...')
    animals = pd.read_csv("./summarised_animals.txt", sep='\t')
    print(animals.columns)
    print(animals['Species'].unique())

    deleted_rows = abundance_plot_utils.delete_rows_for_taxa(
        animals, 'Species', 'Orcinus orca')
    print(deleted_rows)

    # Run the unit tests
    unittest.main()

