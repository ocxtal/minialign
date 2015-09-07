#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>

#define MM_VERSION "0.1.0-r18"

void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

int mm_verbose = 3;
double mm_realtime0;

int main_index(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	int ret = 0, i;
	double t_start;
	liftrlimit();
	if (argc == 1) {
		fprintf(stderr, "Usage: minimap <command> [arguments]\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  index    create minimap index\n");
		return 1;
	}
	mm_realtime0 = t_start = realtime();
	if (strcmp(argv[1], "index") == 0) ret = main_index(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_start, cputime());
	}
	return ret;
}
