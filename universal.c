
/**
 * @file universal.c
 *
 * @brief main function dispatcher
 */
#include <stdio.h>
#include "arch/arch.h"

int sse41_main(int argc, char *argv[]);
int avx2_main(int argc, char *argv[]);


int main(int argc, char *argv[], char *envp[])
{
	#if defined(__x86_64__) && defined(__AVX2__)
	if((arch_cap() & ARCH_CAP_AVX2) != 0) {
		return(avx2_main(argc, argv));
	}
	#endif

	#if defined(__x86_64__) && defined(__SSE4_1__)
	if((arch_cap() & ARCH_CAP_SSE41) != 0) {
		return(sse41_main(argc, argv));
	}
	#endif

	fprintf(stderr, "[E::main] no main function found.\n");
	return(2);
}

/**
 * end of universal.c
 */
