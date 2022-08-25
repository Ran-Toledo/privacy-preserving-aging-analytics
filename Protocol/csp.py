# IMPORTS
import logging
import numpy as np
from phe import paillier
from numpy.linalg import inv

# CONSTANTS
CSP_LEVEL = 8


class CSP:

    def __init__(self, logger):
        self.w = []                     # w~ = masked model
        self.C = []                     # Enc_pk(C)
        self.C_tag = []                 # decrypted C
        self.d = []                     # Enc_pk(d)
        self.d_tag = []                 # decrypted d
        self.t = []                     # t' = masked t vector
        self.__public_key, self.__secret_key = paillier.generate_paillier_keypair()     # public key, private key
        logging.setLoggerClass(logger)

    def get_pk(self):
        return self.__public_key

    def receive_data_from_mle_site_step(self, c: np.array, d: np.array):
        self.C = c
        self.d = d
        logging.log(CSP_LEVEL, "Received C and d.")

    def receive_data_from_mle_time_step(self, t):
        self.t = t
        logging.log(CSP_LEVEL, "Received t.")

    def decrypt(self, num: paillier.EncryptedNumber):
        return self.__secret_key.decrypt(num)

    def decrypt_c_and_d(self):
        self.C_tag = []
        self.d_tag = []

        for enc_row in self.C:
            dec_row = []
            for enc_num in enc_row:
                # if value is encrypted
                if not (isinstance(enc_num, np.integer) or isinstance(enc_num, np.float64)):
                    # Decrypt value and append it
                    dec_row.append(self.decrypt(enc_num))
                else:
                    # append value
                    dec_row.append(enc_num)
            dec_row = np.array(dec_row)
            self.C_tag.append(dec_row)
        self.C_tag = np.array(self.C_tag)

        for enc_num in self.d:
            dec_num = self.decrypt(enc_num)
            self.d_tag.append(dec_num)
        self.d_tag = np.array(self.d_tag)

        logging.log(CSP_LEVEL, "Decrypting C and d complete.")

    def decrypt_t(self):
        # decrypt all values in self.t
        t_tag = [(self.decrypt(enc_num)) for enc_num in self.t]
        self.t = np.array(t_tag)

        logging.log(CSP_LEVEL, "Decrypting t complete.")

    def get_w(self):
        # transpose c_tag
        c_tag_transpose = np.transpose(self.C_tag)

        # compute w~ = masked model
        self.w = np.matmul(np.matmul(inv(np.matmul(c_tag_transpose, self.C_tag)), c_tag_transpose), self.d_tag)

        logging.log(CSP_LEVEL, "Finished computing w.")

    def send_w(self, my_mle):
        logging.log(CSP_LEVEL, "Sending masked w.")
        my_mle.receive_masked_w(self.w)

    def send_t(self, my_mle):
        logging.log(CSP_LEVEL, "Sending masked t.")
        my_mle.receive_masked_t(self.t)
