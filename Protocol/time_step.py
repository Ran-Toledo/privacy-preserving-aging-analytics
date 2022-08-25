# IMPORTS
from csp import CSP
from mle import MLE


def time_step(my_csp: CSP, my_mle: MLE):

    # Step 1
    # MLE uses the output w^* from Phase 2-A to retrieve the rates and initial states
    my_mle.get_rates_and_states()
    my_mle.calc_squared_sum_of_rates()

    # Step 2
    # MLE calculates uses these rates and initial states, and the encrypted datasets to
    # calculate obtain the encryption of the t vector
    my_mle.calc_t_vector()

    # Step 3
    # MLE uses the ciphertexts Enc_pk(t) and private random values in order to obtain
    # encryptions of new values called "masked data"
    my_mle.mask_data_time_step()

    # Step 4
    # These encryptions are then sent to the CSP, which decrypts them.
    my_mle.send_data_to_csp_time_step(my_csp)
    my_csp.decrypt_t()

    # Step 5
    # This computations output is the “masked model”, a vector t', which is then sent back
    # from the CSP to the MLE
    my_csp.send_t(my_mle)

    # Step 6
    # MLE computes the output from t' ( by computing inv(R) * t' )
    my_mle.unmask_t()

    # Step 7
    # This output represents the predicted ages at the end of the current iteration.
    # The ages from the current iteration (self.all_ages) are replaced by these predicted ages.
    my_mle.set_updated_ages()
    my_mle.reset_variables()
