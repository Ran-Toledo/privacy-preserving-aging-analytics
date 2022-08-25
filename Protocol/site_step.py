# IMPORTS
from csp import CSP
from mle import MLE


def site_step(my_csp: CSP, my_mle: MLE):

    # Step 1
    # MLE uses the ciphertexts Enc_pk (A) and Enc_pk (b) and private random values in order
    # to obtain encryptions of new values called “masked data”
    my_mle.get_x_and_y()
    my_mle.mask_data_site_step()

    # Step 2
    # These encryptions are then sent to the CSP, which decrypts and runs a given algorithm
    # on the masked data
    my_mle.send_data_to_csp_site_step(my_csp)
    my_csp.decrypt_c_and_d()
    my_csp.get_w()

    # Step 3
    # This computation's output is the “masked model”, a vector w~, which is then sent back
    # from the CSP to the MLE
    my_csp.send_w(my_mle)

    # Step 4
    # MLE computes the output w^* from w~
    my_mle.unmask_w()
