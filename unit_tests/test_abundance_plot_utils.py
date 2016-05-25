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
        animal_df = load_animal_data()
        self.assertTrue(taxa_name in animal_df[taxa_level].unique())

        # Now remove that taxa and make sure it is gone.
        animal_df = abundance_plot_utils.delete_rows_for_taxa(
        animals, 'Species', 'Orcinus orca')
        self.assertFalse(taxa_name in animal_df[taxa_level].unique())
        print("unique values left: {}".format(animal_df[taxa_level].unique()))

    def test_invalid_taxa_dict(self):
        # !! Doesn't check for all possible modes of invalidity !!
        # load in the animal dict.
        animal_df = pd.read_csv("./summarised_animals.txt", sep='\t')
        # Use a taxa dict that is invalid:
        # Genus	Species
        # Orcinus	Orcinus orca
        invalid_taxa_dict = {'Genus': ['Orcinus'],
                             'Species': ['Orcinus orca']}
        # Check that aggregate_mixed_taxonomy fails, which it should because
        # the values in the taxa dict have overlapping taxa.
        with self.assertRaises(AssertionError):
            # If the code below returns the expected assertion,
            # Python continues on.
            ag = abundance_plot_utils.aggregate_mixed_taxonomy(
            dataframe=animal_df, taxa_dict=invalid_taxa_dict,
            main_dir='../', summarise_other=True)
            # The resulting df has a row of NaN values.
            print(ag)




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

