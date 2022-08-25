# IMPORTS
import utilities as utils
import numpy as np
import logging
from numpy.linalg import inv

# CONSTANTS
MLE_LEVEL = 9
AGES_SILOS_INDEX = 0
MV_SILOS_INDEX = 1


class MLE:

    def __init__(self, pk, logger):
        self.__public_key = pk          # public key
        self.ss_of_rates = None         # squared sum of rates
        self.rates = []                 # rates
        self.initial_states = []        # initial states
        self.all_meth_values = []       # unified encrypted methylation values
        self.__all_ages = []            # unified ages
        self.__original_ages = []       # original encrypted ages
        self.R_matrix = []              # random matrix
        self.r_vec = []                 # random vector
        self.X = []                     # X matrix
        self.y = []                     # y vector
        self.w = []                     # w vector = rates, initial states
        self.C = []                     # C = masked X matrix
        self.d = []                     # d = masked y vector
        self.__t = []                   # t vector
        self.encrypted_data = []        # all separate encrypted datasets
        logging.setLoggerClass(logger)

    def receive_data_from_owners(self, encrypted_data):
        self.encrypted_data.append(encrypted_data)

    def send_data_to_csp_site_step(self, my_csp):
        logging.log(MLE_LEVEL, "Sending C and d.")
        my_csp.receive_data_from_mle_site_step(self.C, self.d)

    def send_data_to_csp_time_step(self, my_csp):
        logging.log(MLE_LEVEL, "Sending t.")
        my_csp.receive_data_from_mle_time_step(self.__t)

    def mask_data_site_step(self):
        rows, cols = np.array(self.X).shape

        # get random matrix and random vector
        self.R_matrix = utils.get_random_matrix(cols)
        self.r_vec = utils.get_random_vector(cols)

        # obtain masked data using random matrix and vector
        self.C = np.matmul(self.X, self.R_matrix)
        self.d = self.y + np.matmul(self.X, self.r_vec)

        logging.log(MLE_LEVEL, "Masking C and d complete.")

    def mask_data_time_step(self):
        rows = len(self.__t)

        # get random matrix
        self.R_matrix = utils.get_random_matrix(rows)

        # obtain masked t vector using random matrix
        self.__t = np.matmul(self.R_matrix, self.__t)

        logging.log(MLE_LEVEL, "Masking t complete.")

    def merge_data(self):

        # Merging datasets
        ages_silos = [dataset[AGES_SILOS_INDEX] for dataset in self.encrypted_data]
        mv_silos = [dataset[MV_SILOS_INDEX] for dataset in self.encrypted_data]

        self.__all_ages = np.concatenate(ages_silos)
        self.__original_ages = self.__all_ages
        self.all_meth_values = np.concatenate(mv_silos, axis=1)

        logging.log(MLE_LEVEL, "Finished merging the datasets.")

    def get_x_and_y(self):

        # Defining X matrix
        n, m = self.all_meth_values.shape
        self.X = np.zeros((m * n, 2 * n), dtype=object)
        row_index = 0
        for j in range(0, n):
            for i in range(0, m):
                self.X[row_index][j] = self.__all_ages[i]
                self.X[row_index][j + n] = 1
                row_index += 1

        # Defining y vector
        for meth_site in self.all_meth_values:
            for gsm_value in meth_site:
                self.y.append(gsm_value)
        self.y = np.array(self.y, dtype=object)

        logging.log(MLE_LEVEL, "Creating X and y complete.")

    def receive_masked_w(self, masked_w):
        logging.log(MLE_LEVEL, "Received masked w.")
        self.w = masked_w

    def receive_masked_t(self, masked_t):
        logging.log(MLE_LEVEL, "Received masked t.")
        self.__t = masked_t

    def unmask_w(self):
        # unmask the model by computing (R*w - r)
        self.w = np.matmul(self.R_matrix, self.w) - self.r_vec

        logging.log(MLE_LEVEL, "Finished unmasking w.")

    def unmask_t(self):
        # unmask the t vector by computing (inv(R)*t)
        self.__t = np.matmul(inv(self.R_matrix), self.__t)

        logging.log(MLE_LEVEL, "Finished unmasking t.")

    def get_rates_and_states(self):
        w_len = len(self.w)
        n = int(w_len / 2)
        self.rates = self.w[:n]
        self.initial_states = self.w[n:]

    def calc_squared_sum_of_rates(self):
        self.ss_of_rates = sum(rate ** 2 for rate in self.rates)

        logging.log(MLE_LEVEL, "Finished calculating squared sum of rates.")

    def calc_t_vector(self):
        self.__t = []
        n, m = self.all_meth_values.shape

        for j in range(0, m):
            t_j = 0
            for i in range(0, n):
                t_j += (self.rates[i] * (self.all_meth_values[i][j] - self.initial_states[i]))
            self.__t.append(t_j / self.ss_of_rates)

        self.__t = np.array(self.__t, dtype=object)

        logging.log(MLE_LEVEL, "Finished calculating t vector.")

    def set_updated_ages(self):
        self.__all_ages = self.__t

    def reset_variables(self):
        self.X = []                 # reset X
        self.y = []                 # reset y
        self.C = []                 # reset C
        self.d = []                 # reset d
        self.__t = []               # reset t
        self.w = []                 # reset w
        self.ss_of_rates = None     # reset squared sum of rates
        self.initial_states = []    # reset initial states
        self.rates = []             # reset rates
        self.R_matrix = []          # reset random matrix
        self.r_vec = []             # reset random vector

    def get_predicted_ages(self):
        return self.__all_ages

    def get_original_ages(self):
        return self.__original_ages
