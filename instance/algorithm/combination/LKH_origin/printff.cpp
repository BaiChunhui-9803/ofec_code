#include <stdio.h>
#include <stdarg.h>
#include "./INCLUDE/LKH.h"
namespace LKH {
	void printff(const char *fmt, ...);

	/*
	 * The printff function prints a message and flushes stdout.
	 */

	void LKHAlg::printff(const char *fmt, ...)
	{
		va_list args;

		va_start(args, fmt);
		vprintf(fmt, args);
		va_end(args);
		fflush(stdout);
	}
}