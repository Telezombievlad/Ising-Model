# Grid Server
# No Copyright. Vladislav Aleinik 2020

import socket
import threading
import socketserver

from enum import Enum

class ClientRequestType(Enum):
	STILL_ALIVE       = 0
	GIVE_ME_TASK      = 1
	HERE_IS_YOUR_DATA = 2

class GridRequest(socketserver.BaseRequestHandler):
	def __init__(self):
		self.type = 'still_alive'

	def __init__(self, temperature, field):
		self.type = "data_request"
		self.temperature = temperature
		self.field = field

	def handle(self):
		



if __name__ == "__main__":

