import numpy as np
import matplotlib.pyplot as plt

from phe import paillier
from numpy.linalg import inv
from scipy import stats
from EpigeneticPacemaker.ExampleData.DataSets import get_example_data
from EpigeneticPacemaker.EpigeneticPacemaker import EpigeneticPacemaker

PRECISION = None
NUMBER_OF_SITES = 5
NUMBER_OF_PEOPLE = 24
MIN_COEFFICIENT = 0.90
MIN_VALUE = 1
MAX_VALUE = 2048


def get_random_matrix(n) -> np.array:
    is_invertible = 1
    r_matrix = np.random.randint(MIN_VALUE, MAX_VALUE, size=(n, n))
    r_inv = inv(r_matrix)

    if np.array_equal(np.matmul(r_matrix, r_inv), np.matmul(r_inv, r_matrix)):
        is_invertible = 0

    while is_invertible == 0:
        is_invertible = 1
        r_matrix = np.random.randint(MIN_VALUE, MAX_VALUE, size=(n, n))

        if np.array_equal(np.matmul(r_matrix, r_inv), np.matmul(r_inv, r_matrix)):
            is_invertible = 0

    return r_matrix


def get_random_vector(n) -> np.array:
    return np.random.randint(MIN_VALUE, MAX_VALUE, size=n)


def pearson_correlation(meth_matrix: np.array, phenotype: np.array) -> np.array:
    # calculate mean for each row and phenotype mean
    matrix_means = np.mean(meth_matrix, axis=1)
    phenotype_mean = np.mean(phenotype)

    # subtract means from observed values
    transformed_matrix = meth_matrix - matrix_means.reshape([-1, 1])
    transformed_phenotype = phenotype - phenotype_mean

    # calculate covariance
    covariance = np.sum(transformed_matrix * transformed_phenotype, axis=1)
    variance_meth = np.sqrt(np.sum(transformed_matrix ** 2, axis=1))
    variance_phenotype = np.sqrt(np.sum(transformed_phenotype ** 2))

    return covariance / (variance_meth * variance_phenotype)


def split_data_to_different_owners(m) -> np.array:
    # retrieve the data
    data_set1, data_set2 = get_example_data()
    samples, cpg_sites, ages, methylation_values = data_set2

    # get the absolute value of the correlation coefficient
    abs_pcc_coefficients = abs(pearson_correlation(methylation_values, ages))

    # return list of site indices with a high absolute correlation coefficient
    training_sites = np.where(abs_pcc_coefficients > MIN_COEFFICIENT)[0]
    training_sites = training_sites[:NUMBER_OF_SITES]
    ages = ages[:NUMBER_OF_PEOPLE]
    methylation_values = methylation_values[training_sites, :NUMBER_OF_PEOPLE]

    methylation_values_split = np.hsplit(methylation_values, m)
    ages = np.array(ages)
    ages_split = np.split(ages, m)

    samples = np.array(samples)
    samples_split = np.split(samples, m)

    data_owners = []
    for i in range(0, m):
        data_owners.append([samples_split[i].tolist(), cpg_sites, ages_split[i], methylation_values_split[i]])

    print("***Splitting data complete.")

    return data_owners


def encrypt_data(data_owner, public_key):
    samples, cpg_sites, ages, methylation_values = data_owner
    cipher_ages = []
    cipher_meth_values = []

    for gsm_age in ages:
        cipher_ages.append(public_key.encrypt(gsm_age, PRECISION))

    print("Encrypting ages complete.")

    for site in methylation_values:
        cipher_site = []
        for gsm_site_value in site:
            cipher_site.append(public_key.encrypt(gsm_site_value, PRECISION))
        cipher_meth_values.append(cipher_site)

    print("Encrypting methylation site rates complete.")

    return [cipher_ages, cipher_meth_values]


def create_graph(x, y):
    plt.scatter(x, y, color="red", s=15)
    plt.xlabel("Chronological ages")
    plt.ylabel("Epigenetic ages")
    plt.gca().xaxis.grid(True)
    plt.gca().yaxis.grid(True)
    res = stats.linregress(x, y)
    plt.plot(np.array(x), res.intercept + res.slope * np.array(x), 'b', label=f"R-squared: {res.rvalue ** 2:.6f}")
    z = np.polyfit(x.flatten(), np.array(y).flatten(), 1)
    # plt.title("y=%.6fx+%.6f"%(z[0],z[1]))
    plt.title("Epigenetic PaceMaker")
    # plt.legend()
    plt.show()
