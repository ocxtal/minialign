
/**
 * @file universal.c
 *
 * @brief main function dispatcher
 */
#include <string.h>

int main_silent(int argc, char *argv[]);
int main_verbose(int argc, char *argv[]);

int main(int argc, char *argv[], char *envp[])
{
	if(argc > 1 && strcmp(argv[1], "verbose") == 0) {
		return(main_verbose(argc - 1, argv + 1));
	} else if(argc > 1 && strcmp(argv[1], "silent") == 0) {
		return(main_silent(argc - 1, argv + 1));
	}
	return(main_silent(argc, argv));
}

/**
 * end of universal.c
 */
