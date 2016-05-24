import unittest

import pandas as pd

class testCompletenessOfSummarisedData(unittest.TestCase):
    def test_animal_data(self):
        animal_df = pd.read_csv("./summarised_animals.txt", sep='\t')
        sum_by_sample = animal_df.groupby(
            ['oxy', 'rep', 'week'])['fraction of reads'].sum()
        self.assertTrue((sum_by_sample > 0.999).all())
        self.assertTrue((sum_by_sample < 1.001).all())


class testAbundannceSummary(unittest.TestCase):
    def test_summary_with_all_taxonomy_remaining(self):
        summary_df = \
            pd.read_csv("../results/reduced_data--all_taxonomy_remains.csv")
        sum_by_sample = summary_df.groupby(
            ['oxy', 'rep', 'week'])['fraction of reads'].sum()
        self.assertTrue((sum_by_sample > 0.999).all())
        self.assertTrue((sum_by_sample < 1.001).all())


if __name__ == '__main__':

    animal_df = pd.read_csv("./summarised_animals.txt", sep='\t')
    print(animal_df.head())
    sums = animal_df.groupby(
            ['oxy', 'rep', 'week'])['fraction of reads'].sum()
    # make sure all the sums are 1:
    print((sums > 0.99).all())
    print((sums < 1.01).all())
    print(sums)


    # Run the Unit Tests
    unittest.main()
