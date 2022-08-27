# IMPORTS
import numpy as np
import utilities as utils
from EpigeneticPacemaker.ExampleData.DataSets import get_example_data


# CONSTANTS
PRECISION = None
MIN_COEFFICIENT = 0.90
NUMBER_OF_SITES = 5
NUMBER_OF_PEOPLE = 24


def split_data_to_different_owners(m: int) -> np.array:

    # Retrieve the data
    data_set1, data_set2 = get_example_data()
    samples, cpg_sites, ages, methylation_values = data_set2

    assert (len(ages) % m == 0), "Cannot split data into separate owners"

    # Get the absolute value of the correlation coefficient
    abs_pcc_coefficients = abs(utils.pearson_correlation(methylation_values, ages))

    # Return list of site indices with a high absolute correlation coefficient
    training_sites = np.where(abs_pcc_coefficients > MIN_COEFFICIENT)[0]
    training_sites = training_sites[:NUMBER_OF_SITES]
    ages = ages[:NUMBER_OF_PEOPLE]
    methylation_values = methylation_values[training_sites, :NUMBER_OF_PEOPLE]

    # Split data
    methylation_values_split = np.hsplit(methylation_values, m)
    ages = np.array(ages)
    ages_split = np.split(ages, m)
    samples = np.array(samples)
    samples_split = np.split(samples, m)

    # All separated datasets
    data_owners = ([samples_split[i].tolist(), cpg_sites, ages_split[i], methylation_values_split[i]] for i in range(m))

    return data_owners


def encrypt_data(data_owner, public_key):

    # Unpack data owner's dataset
    samples, cpg_sites, ages, methylation_values = data_owner

    # Encrypt ages
    cipher_ages = [public_key.encrypt(age, PRECISION) for age in ages]

    # Encrypt methylation values
    cipher_meth_values = []
    for site in methylation_values:
        encrypted_site = [(public_key.encrypt(gsm_site_value, PRECISION)) for gsm_site_value in site]
        cipher_meth_values.append(encrypted_site)

    return [cipher_ages, cipher_meth_values]
