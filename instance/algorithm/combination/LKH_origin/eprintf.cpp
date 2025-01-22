#include "./INCLUDE/LKH.h"
#include <stdarg.h>
#include "../../../../core/exception.h"


/* 
 * The eprintf function prints an error message and exits.
 */

void LKH::LKHAlg::eprintf(const char *fmt, ...)
{
    std::string lkherror = "lkh error\t" + m_filename + "\t" + fmt + "\n";
    throw ofec::Exception(lkherror.c_str());

   /* va_list args;

    if (LastLine && *LastLine)
        fprintf(stderr, "\n%s\n", LastLine);
    fprintf(stderr, "\n*** Error ***\n");
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);*/
}
