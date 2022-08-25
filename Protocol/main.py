# IMPORTS
import utilities as utils
import logging
import time
from csp import CSP
from mle import MLE
from site_step import site_step
from time_step import time_step

# CONSTANTS
NUM_OF_OWNERS = 8
NUM_OF_ITERATIONS = 4
CSP_LEVEL = 8
MLE_LEVEL = 9

# ABSTRACT
"""
Phase1: 	
    1.      CSP generates the key pair (sk,pk), stores sk and makes pk public
    2.      Each data owner encrypts their dataset using pk and sends it to the MLE
    3.      MLE uses these encrypted datasets and the homomorphic property of the encryption scheme
            in order to obtain encryptions of A and b (coefficient matrix and vector)
            A = X   ,   b = y   (X, y as defined in the EPM paper)
            
Phase2-A:
    1.     	MLE uses the ciphertexts Enc_pk (A) and Enc_pk (b) and private random values in order
            to obtain encryptions of new values called “masked data”
    2.      These encryptions are then sent to the CSP, which decrypts and runs a given algorithm
            on the masked data
    3.      This computation's output is the “masked model”, a vector w~, which is then sent back
            from the CSP to the MLE
    4.      MLE computes the output w^* from w~

Phase2-B:
    1.     	MLE uses the output w^* from Phase 2-A to retrieve the rates and initial states
    2.      MLE calculates uses these rates and initial states, and the encrypted datasets to 
            calculate obtain the encryption of the t vector
    3.      MLE uses the ciphertexts Enc_pk(t) and private random values in order to obtain
            encryptions of new values called "masked data"
    4.      These encryptions are then sent to the CSP, which decrypts them.
    5.      This computations output is the “masked model”, a vector t', which is then sent back
            from the CSP to the MLE
    6.      MLE computes the output from t' ( by computing inv(R) * t' )
    7.      This output represents the predicted ages at the end of the current iteration.
            The ages from the current iteration (self.all_ages) are replaced by these predicted ages.
"""

# MAIN
if __name__ == '__main__':

    # Config logging
    logger = utils.config_logging(logging.getLoggerClass())
    logging.setLoggerClass(logger)

    # Init CSP
    MyCSP = CSP(logger)
    public_key = MyCSP.get_pk()

    # Init MLE
    MyMLE = MLE(public_key, logger)

    # Split dataset obtained from the Epigenetic PaceMaker module
    data_owners = utils.split_data_to_different_owners(NUM_OF_OWNERS)

    # Data Encryption
    for data_owner in data_owners:

        # Each data owner encrypts their dataset
        encrypted_data = utils.encrypt_data(data_owner, public_key)

        # The encrypted dataset is sent to the MLE
        MyMLE.receive_data_from_owners(encrypted_data)

    logging.info("Data Encryption complete.")

    # MLE merges the encrypted datasets it has received from the data owners
    MyMLE.merge_data()

    logging.info("Merging data complete.")

    logging.info("Starting protocol.")

    # Log start time of protocol
    start_time = time.time()

    # Protocol
    for iteration in range(NUM_OF_ITERATIONS):
        logging.info("BEGINNING ITERATION " + str(iteration+1))

        # Perform site step
        site_step(MyCSP, MyMLE)
        logging.info("SITE STEP COMPLETE")

        # Perform time step
        time_step(MyCSP, MyMLE)
        logging.info("TIME STEP COMPLETE")

    logging.info("Completed protocol.")

    end_time = time.time()

    hours, minutes, seconds = utils.timer(start_time, end_time)

    logging.info("Runtime: " + ("{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)))

    # Print Results
    original_ages = [MyCSP.decrypt(age) for age in MyMLE.get_original_ages()]
    utils.create_graph(original_ages, MyMLE.get_predicted_ages())

