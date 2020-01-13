// No Copyright. Vladislav Aleinik 2019

#include <sys/socket.h>
#include <netinet/ip.h>
#include <sys/types.h>
#include <netdb.h>

int main()
{
	// ========================================================================
	// Aquire And Bind Connection Socket                                                
	// ========================================================================
	
	int server_sock = socket(AF_INET, SOCK_STREAM, 0);
	if (server_sock == -1)
	{
		fprintf(stderr, "[COMPUTING-GRID] Unable to get socket\n");
		exit(EXIT_FAILURE);
	}

	struct addrinfo hints;
	struct addrinfo* addr;

	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;

	// ========================================================================
	// Server Accept Loop                                                
	// ========================================================================

	while (true)
	{

	}

	return EXIT_SUCCESS;
}