// No Copyright. Vladislav Aleinik 2019

int main()
{
	// ========================================================================
	// Aquire Connection Socket                                                
	// ========================================================================
	
	int conn_sock = socket(AF_INET, SOCK_STREAM, 0);
	if (conn_sock == -1)
	{
		fprintf(stderr, "[COMPUTING-GRID] Unable to get socket\n");
		exit(EXIT_FAILURE);
	}

	return EXIT_SUCCESS;
}