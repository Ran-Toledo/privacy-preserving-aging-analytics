import numpy as np
from phe import paillier
from numpy.linalg import inv

'''
import socket

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((socket.gethostname(), 6060))
s.listen(5)

while True:
    mleSocket, address = s.accept()
    print("Connection with MLE established from address:", address)
    mleSocket.send(bytes("This is a message from the CSP", "utf-8"))
'''
'''
import socket
import threading

HEADER = 64
PORT = 6060
SERVER = socket.gethostname()
ADDR = (SERVER, PORT)
FORMAT = 'utf-8'
DISCONNECT_MESSAGE = "!DISCONNECT"


MyCSP = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
MyCSP.bind(ADDR)


def handle_client(conn, addr):
    print(f"[NEW CONNECTION] {addr} connected.")

    connected = True
    while connected:
        msg_length = conn.recv(HEADER).decode(FORMAT)
        if msg_length:
            msg_length = int(msg_length)
            msg = conn.recv(msg_length).decode(FORMAT)
            if msg == DISCONNECT_MESSAGE:
                connected = False

            print(f"[MLE] {msg}")
            conn.send("Msg received".encode(FORMAT))

    conn.close()


def start():
    MyCSP.listen()
    print(f"[LISTENING] Server is listening on {SERVER}")
    while True:
        conn, addr = MyCSP.accept()
        conn.send("Connection established.".encode(FORMAT))
        thread = threading.Thread(target=handle_client, args=(conn, addr))
        thread.start()
        print(f"[ACTIVE CONNECTIONS] {threading.activeCount() - 1}")


print("[STARTING] server is starting...")
start()
'''


class CSP:

    def __init__(self):
        self.w = []
        self.C = []
        self.C_tag = []
        self.d = []
        self.d_tag = []
        self.t = []
        self.t_tag = []
        self.__public_key, self.__secret_key = paillier.generate_paillier_keypair()

    def get_pk(self):
        return self.__public_key

    def receive_data_from_mle(self, C, D):
        self.C = C
        self.d = D
        print("CSP:\t\t\tReceived C and d.")

    def receive_data_from_mle_time_step(self, t):
        self.t = t
        print("CSP:\t\t\tReceived t.")

    def decrypt(self, num):
        return self.__secret_key.decrypt(num)

    def decrypt_c_and_d(self):
        self.C_tag = []
        self.d_tag = []

        print("CSP:\t\t\tDecrypting C and d....")

        for enc_row in self.C:
            dec_row = []
            for enc_num in enc_row:
                if not (isinstance(enc_num, np.integer) or isinstance(enc_num, np.float64)):
                    dec_num = self.decrypt(enc_num)
                    dec_row.append(dec_num)
                else:
                    dec_row.append(enc_num)
            dec_row = np.array(dec_row)
            self.C_tag.append(dec_row)
        self.C_tag = np.array(self.C_tag)

        for enc_num in self.d:
            dec_num = self.decrypt(enc_num)
            self.d_tag.append(dec_num)
        self.d_tag = np.array(self.d_tag)

        print("CSP:\t\t\tDecrypting C and d complete.")

    def decrypt_t(self):
        self.t_tag = []

        print("CSP:\t\t\tDecrypting t....")

        for enc_num in self.t:
            dec_num = self.decrypt(enc_num)
            self.t_tag.append(dec_num)
        self.t_tag = np.array(self.t_tag)

        print("CSP:\t\t\tDecrypting t complete.")

    def get_w(self):
        print("CSP:\t\t\tComputing w....")

        c_tag_transpose = np.transpose(self.C_tag)
        self.w = np.matmul(np.matmul(inv(np.matmul(c_tag_transpose, self.C_tag)), c_tag_transpose), self.d_tag)

        print("CSP:\t\t\tFinished computing w.")

    def send_w(self, my_mle):
        print("CSP:\t\t\tSending masked w.")
        my_mle.receive_masked_w(self.w)

    def send_t(self, my_mle):
        print("CSP:\t\t\tSending masked t.")
        my_mle.receive_masked_t(self.t_tag)
