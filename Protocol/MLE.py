'''
import socket

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((socket.gethostname(), 6060))

message = s.recv(2048)

print("Message received:", message)
'''

import socket

HEADER = 64
PORT = 6060
FORMAT = 'utf-8'
DISCONNECT_MESSAGE = "!DISCONNECT"
SERVER = socket.gethostname()
ADDR = (SERVER, PORT)

MLEserver = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
MLEserver.connect(ADDR)


def send(msg):
    message = msg.encode(FORMAT)
    msg_length = len(message)
    send_length = str(msg_length).encode(FORMAT)
    send_length += b' ' * (HEADER - len(send_length))
    MLEserver.send(send_length)
    MLEserver.send(message)
    msg = MLEserver.recv(2048).decode(FORMAT)
    print(f"[CSP] {msg}")


while True:
    msg = MLEserver.recv(2048).decode(FORMAT)
    print(f"[CSP] {msg}")
