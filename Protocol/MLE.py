import utilities as utils
import numpy as np
from numpy.linalg import inv

'''
import socket

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((socket.gethostname(), 6060))

message = s.recv(2048)

print("Message received:", message)
'''
'''
import socket

HEADER = 64
PORT = 6060
FORMAT = 'utf-8'
DISCONNECT_MESSAGE = "!DISCONNECT"
SERVER = socket.gethostname()
ADDR = (SERVER, PORT)

MyMLE = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
MyMLE.connect(ADDR)


def send(msg):
    message = msg.encode(FORMAT)
    msg_length = len(message)
    send_length = str(msg_length).encode(FORMAT)
    send_length += b' ' * (HEADER - len(send_length))
    MyMLE.send(send_length)
    MyMLE.send(message)
    msg = MyMLE.recv(2048).decode(FORMAT)
    print(f"[CSP] {msg}")


while True:
    msg = MyMLE.recv(2048).decode(FORMAT)
    print(f"[CSP] {msg}")
'''


class MLE:

    def __init__(self, pk):
        self.__public_key = pk
        self.ss_of_rates = None
        self.rates = []
        self.initial_states = []
        self.w = []
        self.all_meth_values = []
        self.all_ages = []
        self.original_ages = []
        self.R_matrix = []
        self.r_vec = []
        self.X = []
        self.y = []
        self.C = []
        self.d = []
        self.t = []
        self.encrypted_data = []

    def receive_data_from_owners(self, encrypted_data):
        self.encrypted_data.append(encrypted_data)

    def send_data_to_csp_site_step(self, MyCSP):
        print("MLE:\t\t\tSending C and d.")
        MyCSP.receive_data_from_mle(self.C, self.d)

    def send_data_to_csp_time_step(self, MyCSP):
        print("MLE:\t\t\tSending t.")
        MyCSP.receive_data_from_mle_time_step(self.t)

    def mask_data_site_step(self):
        rows, cols = self.X.shape

        self.R_matrix = utils.get_random_matrix(cols)
        self.r_vec = utils.get_random_vector(cols)

        self.C = np.matmul(self.X, self.R_matrix)
        self.d = self.y + np.matmul(self.X, self.r_vec)

        print("MLE:\t\t\tMasking C and d complete.")

    def mask_data_time_step(self):
        rows = len(self.t)

        self.R_matrix = utils.get_random_matrix(rows)

        self.t = np.matmul(self.R_matrix, self.t)

        print("MLE:\t\t\tMasking t complete.")

    def merge_data(self):
        print("MLE:\t\t\tMerging the encrypted data sets received from owners.")

        # Merging the data
        ages_silos = []
        mv_silos = []
        for data in self.encrypted_data:
            ages_silos.append(data[0])
            mv_silos.append(data[1])

        self.all_ages = np.concatenate(ages_silos)
        self.original_ages = self.all_ages
        self.all_meth_values = np.concatenate(mv_silos, axis=1)

        print("MLE:\t\t\tFinished merging the datasets.")

    def get_x_and_y(self):

        # Defining X matrix
        n, m = self.all_meth_values.shape
        self.X = np.zeros((m * n, 2 * n), dtype=object)
        row_index = 0
        for j in range(0, n):
            for i in range(0, m):
                self.X[row_index][j] = self.all_ages[i]
                self.X[row_index][j + n] = 1
                row_index += 1

        # Defining y vector
        for meth_site in self.all_meth_values:
            for gsm_value in meth_site:
                self.y.append(gsm_value)
        self.y = np.array(self.y, dtype=object)

        print("MLE:\t\t\tCreating X and y complete.")

    def receive_masked_w(self, masked_w):
        print("MLE:\t\t\tReceived masked w.")
        self.w = masked_w

    def receive_masked_t(self, masked_t):
        print("MLE:\t\t\tReceived masked t.")
        self.t = masked_t

    def unmask_w(self):
        print("MLE:\t\t\tUnmasking w....")

        self.w = np.matmul(self.R_matrix, self.w) - self.r_vec

        print("MLE:\t\t\tFinished unmasking w.")

    def unmask_t(self):
        print("MLE:\t\t\tUnmasking t....")

        self.t = np.matmul(inv(self.R_matrix), self.t)

        print("MLE:\t\t\tFinished unmasking t.")

    def get_rates_and_states(self):
        double_n = self.w.shape
        n = int(double_n[0] / 2)
        self.rates = self.w[:n]
        self.initial_states = self.w[n:]

    def calc_squared_sum_of_rates(self):
        print("MLE:\t\t\tCalculating squared sum of rates....")

        self.ss_of_rates = sum(rate ** 2 for rate in self.rates)

        print("MLE:\t\t\tFinished calculating squared sum of rates.")

    def calc_tj_vector(self):
        print("MLE:\t\t\tCalculating t_j vector....")
        self.t = []
        n, m = self.all_meth_values.shape

        for j in range(0, m):
            t_j = 0
            for i in range(0, n):
                t_j += (self.rates[i] * (self.all_meth_values[i][j] - self.initial_states[i]))
            self.t.append(t_j / self.ss_of_rates)

        self.t = np.array(self.t, dtype=object)

        print("MLE:\t\t\tFinished calculating t_j vector.")

    def set_updated_ages(self):
        self.all_ages = self.t

    def reset_variables(self):
        self.X = []
        self.y = []
        self.C = []
        self.d = []
        self.t = []
        self.w = []
        self.ss_of_rates = None
        self.initial_states = []
        self.rates = []
        self.R_matrix = []
        self.r_vec = []

    def get_predicted_ages(self):
        return self.all_ages
