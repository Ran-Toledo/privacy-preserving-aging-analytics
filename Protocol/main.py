# IMPORTS
import utilities as utils
import data_handler
import logging
import time
import epigenetic_pacemaker
from csp import CSP
from mle import MLE
from site_step import site_step
from time_step import time_step


# CONSTANTS
NUM_OF_OWNERS = 8
NUM_OF_ITERATIONS = 4


# ABSTRACT
"""
Phase1: 	
    1.      CSP generates the key pair (sk,pk), stores sk and makes pk public
    2.      Each data owner encrypts their dataset using pk and sends it to the MLE
            
Phase2-A:
    1.      MLE uses these encrypted datasets and the homomorphic property of the encryption scheme
            in order to obtain encryptions of X and y (as defined in the EPM paper)
    2.     	MLE uses the ciphertexts Enc_pk (X) and Enc_pk (y) and private random values in order
            to obtain encryptions of new values called “masked data”
    3.      These encryptions are then sent to the CSP, which decrypts and runs a given algorithm
            on the masked data
    4.      This computation's output is the “masked model”, a vector w~, which is then sent back
            from the CSP to the MLE
    5.      MLE computes the output w^* from w~

Phase2-B:
    1.     	MLE uses the output w^* from Phase 2-A to retrieve the rates and initial states
    2.      MLE calculates uses these rates and initial states, and the encrypted datasets to 
            calculate obtain the encryption of the t vector
    3.      MLE uses the ciphertexts Enc_pk(t) and private random values in order to obtain
            encryptions of new values called "masked data"
    4.      This encrypted and masked vector is then sent to the CSP, which decrypts it.
    5.      The output is the “masked model”, a vector t', which is then sent back
            from the CSP to the MLE
    6.      MLE computes the output from t' ( by computing inv(R) * t' )
    7.      This output represents the predicted ages at the end of the current iteration.
            The ages from the current iteration (self.all_ages) are replaced by these predicted ages.
"""

# MAIN
if __name__ == '__main__':

    # Config logging
    logger = utils.config_logging()
    logging.setLoggerClass(logger)

    # Init CSP
    my_csp = CSP(logger)
    public_key = my_csp.get_pk()

    # Init MLE
    my_mle = MLE(public_key, logger)

    # Split dataset obtained from the Epigenetic PaceMaker module
    data_owners = data_handler.split_data_to_different_owners(NUM_OF_OWNERS)

    logging.info("Starting protocol.")

    # Log start time of protocol
    start_time = time.time()

    # Data Encryption
    for data_owner in data_owners:

        # Each data owner encrypts their dataset
        encrypted_data = data_handler.encrypt_data(data_owner, public_key)

        # The encrypted dataset is sent to the MLE
        my_mle.receive_data_from_owners(encrypted_data)

    logging.info("Data Encryption complete.")

    # MLE merges the encrypted datasets it has received from the data owners
    my_mle.merge_data()

    # Protocol
    for iteration in range(NUM_OF_ITERATIONS):
        logging.info("STARTING ITERATION " + str(iteration+1))

        # Perform site step
        site_step(my_csp, my_mle)
        logging.info("SITE STEP COMPLETE")

        # Perform time step
        time_step(my_csp, my_mle)
        logging.info("TIME STEP COMPLETE")

    logging.info("Completed protocol.")

    end_time = time.time()

    hours, minutes, seconds = utils.timer(start_time, end_time)

    logging.info("Runtime: " + ("{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)))

    # Print Results
    original_ages = [my_csp.decrypt(age) for age in my_mle.get_original_ages()]
    utils.create_graph(original_ages, my_mle.get_predicted_ages())

    # Run algorithm on unencrypted dataset
    # predicted_ages = epigenetic_pacemaker.run_unencrypted_epm()
    # mean_error = utils.mean_average_error(predicted_ages, my_mle.get_predicted_ages())

    # logging.info("Mean average error: " + str(mean_error))
