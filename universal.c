
/**
 * @file universal.c
 *
 * @brief main function dispatcher
 */
#define ARCH_CAP

#include <stdio.h>
#include "arch/arch.h"

#if defined(__x86_64__)
	int sse41_main(int argc, char *argv[]);
	int avx2_main(int argc, char *argv[]);
#elif defined(AARCH64)
#elif defined(PPC64)
#endif

int main(int argc, char *argv[], char *envp[])
{
	#if defined(__x86_64__)
	if((arch_cap() & ARCH_CAP_AVX2) != 0) {
		return(avx2_main(argc, argv));
	}
	if((arch_cap() & ARCH_CAP_SSE41) != 0) {
		return(sse41_main(argc, argv));
	}
	#elif defined(AARCH64)
	#elif defined(PPC64)
	#endif

	fprintf(stderr, "[E::main] no main function found.\n");
	return(2);
}

/**
 * end of universal.c
 */
