
/**
 * @file log.h
 *
 * @brief log handler
 */
#ifndef _LOG_H_INCLUDED
#define _LOG_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

/**
 * color outputs
 */
#define RED(x)			"\x1b[31m" \
						x \
						"\x1b[39m"
#define GREEN(x)		"\x1b[32m" \
						x \
						"\x1b[39m"
#define YELLOW(x)		"\x1b[33m" \
						x \
						"\x1b[39m"
#define BLUE(x)			"\x1b[34m" \
						x \
						"\x1b[39m"
#define MAGENTA(x)		"\x1b[35m" \
						x \
						"\x1b[39m"
#define CYAN(x)			"\x1b[36m" \
						x \
						"\x1b[39m"
#define WHITE(x)		"\x1b[37m" \
						x \
						"\x1b[39m"

/**
 * @macro DEBUG
 */
// #define DEBUG 			( 1 )

/**
 * @macro dbprintf
 */
#ifdef DEBUG

	#define dbprintf(fmt, ...) { \
		fprintf(stderr, fmt, __VA_ARGS__); \
	}

#else

	#define dbprintf(fmt, ...) {}

#endif

/**
 * @macro debug
 */
#ifdef DEBUG

	#define debug(...) { \
		debug_impl(__VA_ARGS__, ""); \
	}
	#define debug_impl(fmt, ...) { \
		dbprintf("[%s: %s(%d)] " fmt "%s\n", __FILE__, __func__, __LINE__, __VA_ARGS__); \
	}

#else

	#define debug(...) {}

#endif

/**
 * @macro print_lane
 */
#ifdef DEBUG

	#define print_lane(p1, p2) { \
		cell_t *p = p1, *t = p2; \
		char *str = NULL; \
		int len = 0, size = 256; \
		str = malloc(size); \
		len += sprintf(str+len, "["); \
		while(p != t) { \
			if(*--t <= CELL_MIN) { \
				len += sprintf(str+len, "-oo,"); \
			} else if(*t >= CELL_MAX) { \
				len += sprintf(str+len, "oo,"); \
			} else { \
				len += sprintf(str+len, "%d,", *t); \
			} \
			if(len > (size - 20)) { \
				size *= 2; \
				str = realloc(str, size); \
			} \
		} \
		str[len == 1 ? 1 : len-1] = ']'; \
		str[len == 1 ? 2 : len] = '\0'; \
		debug("lane(%s)", str); \
		free(str); \
	}

#else

	#define print_lane(p1, p2) {}

#endif

/**
 * @macro logprintf
 */
#define logprintf(fmt, ...) { \
	fprintf(stderr, fmt, __VA_ARGS__); \
}

/**
 * @macro log
 */
#define log(...) { \
	log_impl(__VA_ARGS__, ""); \
}
#define log_impl(fmt, ...) { \
	logprintf("[%s] " fmt "%s\n", __func__, __VA_ARGS__); \
}
#define log_nr(...) { \
	log_nr_impl(__VA_ARGS__, ""); \
}
#define log_nr_impl(fmt, ...) { \
	logprintf("[%s] " fmt "%s", __func__, __VA_ARGS__); \
}

/**
 * @macro log_error
 */
#define log_error(...) { \
	log_error_impl(__VA_ARGS__, ""); \
}
#define log_error_impl(fmt, ...) { \
	logprintf("[%s]  ERROR: " fmt "%s\n", __func__, __VA_ARGS__); \
}

/**
 * @macro log_error_abort
 */
#define log_error_abort(...) { \
	log_error_impl(__VA_ARGS__, ""); \
	exit(1); \
}

/**
 * @macro msg
 */
#define msg(...) { \
	msg_impl(__VA_ARGS__, ""); \
}
#define msg_impl(fmt, ...) { \
	logprintf(fmt "%s\n", __VA_ARGS__); \
}

/**
 * @macro dump
 */
#ifndef dump

#if DEBUG
	/* compatible with dump in unittest.h */
	#define ut_dump(ptr, len) ({ \
		uint64_t size = (((len) + 15) / 16 + 1) * \
			(strlen("0x0123456789abcdef:") + 16 * strlen(" 00a") + strlen("  \n+ margin")) \
			+ strlen(#ptr) + strlen("\n`' len: 100000000"); \
		uint8_t *_ptr = (uint8_t *)(ptr); \
		char *_str = alloca(size); \
		char *_s = _str; \
		/* make header */ \
		_s += sprintf(_s, "\n`%s' len: %" PRId64 "\n", #ptr, (int64_t)len); \
		_s += sprintf(_s, "                   "); \
		for(int64_t i = 0; i < 16; i++) { \
			_s += sprintf(_s, " %02x", (uint8_t)i); \
		} \
		_s += sprintf(_s, "\n"); \
		for(int64_t i = 0; i < ((len) + 15) / 16; i++) { \
			_s += sprintf(_s, "0x%016" PRIx64 ":", (uint64_t)_ptr); \
			for(int64_t j = 0; j < 16; j++) { \
				_s += sprintf(_s, " %02x", (uint8_t)_ptr[j]); \
			} \
			_s += sprintf(_s, "  "); \
			for(int64_t j = 0; j < 16; j++) { \
				_s += sprintf(_s, "%c", isprint(_ptr[j]) ? _ptr[j] : ' '); \
			} \
			_s += sprintf(_s, "\n"); \
			_ptr += 16; \
		} \
		(char const *)_str; \
	})
#else

	#define dump(ptr, len) ;

#endif
#endif /* #ifndef dump */

#endif /* #ifndef _LOG_H_INCLUDED */
/**
 * end of log.h
 */
