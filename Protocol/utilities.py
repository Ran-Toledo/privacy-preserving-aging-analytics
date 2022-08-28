# IMPORTS
import sys
import time
import logging
import numpy as np
import matplotlib.pyplot as plt
from logging import Logger
from numpy.linalg import inv


# CONSTANTS
MIN_VALUE = 1
MAX_VALUE = 2048
CSP_LEVEL = 8
MLE_LEVEL = 9


def timer(start: time, end: time) -> (int, int, int):
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    return hours, minutes, seconds


def config_logging() -> [Logger]:
    logging.addLevelName(CSP_LEVEL, "CSP")
    logging.addLevelName(MLE_LEVEL, "MLE")
    logging.basicConfig(
        level=CSP_LEVEL,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.basicConfig(
        level=MLE_LEVEL,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

    return logging.getLoggerClass()


def get_random_matrix(n: int) -> np.array:

    # invertible matrix indicator
    is_invertible = 1

    # sample random matrix and compute its inverse
    r_matrix = np.random.randint(MIN_VALUE, MAX_VALUE, size=(n, n))
    r_inv = inv(r_matrix)

    # check if array is invertible
    if np.array_equal(np.matmul(r_matrix, r_inv), np.matmul(r_inv, r_matrix)):
        is_invertible = 0

    # continue sampling random matrices until invertible
    while is_invertible == 0:
        is_invertible = 1
        r_matrix = np.random.randint(MIN_VALUE, MAX_VALUE, size=(n, n))

        if np.array_equal(np.matmul(r_matrix, r_inv), np.matmul(r_inv, r_matrix)):
            is_invertible = 0

    return r_matrix


def get_random_vector(n: int) -> np.array:

    # sample random vector
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


def create_graph(x: np.array, y: np.array, sites: int, individuals: int):
    logging.disable(logging.DEBUG)
    plt.title("Epigenetic PaceMaker\nSites = {}, Individuals = {}".format(str(sites), str(individuals)))
    plt.xlabel("Chronological ages")
    plt.ylabel("Epigenetic ages")
    plt.xticks([tick for tick in range(-20, 100, 5)])
    plt.yticks([tick for tick in range(-20, 100, 5)])
    plt.gca().xaxis.grid(True)
    plt.gca().yaxis.grid(True)
    plt.scatter(x, y, color="blue", s=5)
    plt.show()
