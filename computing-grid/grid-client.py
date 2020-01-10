# Grid Client
# No Copyright. Vladislav Aleinik 2020  

import socket
import socketserver
import sys

from enum import Enum

class ServerRequestType(Enum):
	STILL_ALIVE       = 0
	HERE_IS_YOUR_TASK = 1
	DATA_DONE         = 2




if __name__ == "__main__":

	