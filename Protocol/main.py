import CSP
import MLE
import utilities as utils

NUM_OF_OWNERS = 8
NUM_OF_ITERATIONS = 4


def site_step(MyCSP, MyMLE):
    """
    Phase2A:
        1.     	MLE uses the ciphertexts Enc_pk (A) and Enc_pk (b) and private random values in order
                to obtain encryptions of new values called “masked data”
        2.      These encryptions are then sent to the CSP, which decrypts and runs a given algorithm
                on the masked data
        3.      This computations output is the “masked model”, a vector w~, which is then sent back
                from the CSP to the MLE
        4.      MLE computes the output w^* from w~

    """
    # Step 1
    MyMLE.get_x_and_y()
    MyMLE.mask_data_site_step()

    # Step 2
    MyMLE.send_data_to_csp_site_step(MyCSP)
    MyCSP.decrypt_c_and_d()
    MyCSP.get_w()

    # Step 3
    MyCSP.send_w(MyMLE)

    # Step 4
    MyMLE.unmask_w()


def time_step(MyCSP, MyMLE):
    """
    Phase2B:
        1.     	MLE: Calculates sum of r_i squared (the same for each t_j)
        2.      MLE: For each t_j , calculates the sum of (r_i * (s_i,j - s0_i) )
        3.      MLE: Uses private randon values in order to ebtain encryptions of new values
                called "masked data" (A=I , b = vector with values of all t_j)
        4.      These encryptions are then sent to the CSP, which decrypts and runs a given algorithm
                on the masked data (w = A^-1*b = A*b = b)
        5.      This computations output is the “masked model”, a vector w~, which is then sent back
                from the CSP to the MLE
        6.      MLE computes the output w^* from w~

    """
    # Step 1
    MyMLE.get_rates_and_states()
    MyMLE.calc_squared_sum_of_rates()

    # Step 2
    MyMLE.calc_tj_vector()

    # Step 3
    MyMLE.mask_data_time_step()

    # Step 4
    MyMLE.send_data_to_csp_time_step(MyCSP)
    MyCSP.decrypt_t()

    # Step 5
    MyCSP.send_t(MyMLE)

    # Step 6
    MyMLE.unmask_t()
    MyMLE.set_updated_ages()
    MyMLE.reset_variables()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    """
    Phase1: 	
        1.      CSP generates the key pair (sk,pk), stores sk and makes pk public
        2.      Each data owner DO_i sends MLE ciphertexts computed using pk and the values in D_i
        3.      MLE uses these ciphertexts and the homomorphic property of the encryption scheme
                in order to obtain encryptions of A and b (coefficient matrix and vector)
                A = X   ,   b = y
    """
    # Step 1
    MyCSP = CSP.CSP()
    public_key = MyCSP.get_pk()

    # Step 2
    MyMLE = MLE.MLE(public_key)
    data_owners = utils.split_data_to_different_owners(NUM_OF_OWNERS)

    for data_owner in data_owners:
        # Each data owner encrypts their dataset
        encrypted_data = utils.encrypt_data(data_owner, public_key)

        # The encrypted dataset is sent to the MLE
        MyMLE.receive_data_from_owners(encrypted_data)

    print("*********************************************")
    print("Data Encryption complete. Beginning protocol.")
    print("*********************************************")

    # Step 3
    MyMLE.merge_data()

    for i in range(0, NUM_OF_ITERATIONS):
        site_step(MyCSP, MyMLE)
        time_step(MyCSP, MyMLE)
