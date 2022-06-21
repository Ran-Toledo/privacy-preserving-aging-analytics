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
