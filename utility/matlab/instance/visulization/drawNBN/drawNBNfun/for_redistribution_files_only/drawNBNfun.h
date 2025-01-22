//
// MATLAB Compiler: 8.2 (R2021a)
// Date: Tue May 30 19:12:58 2023
// Arguments:
// "-B""macro_default""-W""cpplib:drawNBNfun,all,version=1.0""-T""link:lib""-d""
// E:\Diao_Yiya\code\OFEC\utility\matlab\instance\drawNBN\drawNBNfun\for_testing
// ""-v""\\172.24.242.8\share\Student\2018\YiyaDiao\nbn_code\matlab_visual\drawN
// BNfun.m"
//

#ifndef drawNBNfun_h
#define drawNBNfun_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#include "mclcppclass.h"
#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_drawNBNfun_C_API 
#define LIB_drawNBNfun_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_drawNBNfun_C_API 
bool MW_CALL_CONV drawNBNfunInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_drawNBNfun_C_API 
bool MW_CALL_CONV drawNBNfunInitialize(void);

extern LIB_drawNBNfun_C_API 
void MW_CALL_CONV drawNBNfunTerminate(void);

extern LIB_drawNBNfun_C_API 
void MW_CALL_CONV drawNBNfunPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_drawNBNfun_C_API 
bool MW_CALL_CONV mlxDrawNBNfun(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_drawNBNfun
#define PUBLIC_drawNBNfun_CPP_API __declspec(dllexport)
#else
#define PUBLIC_drawNBNfun_CPP_API __declspec(dllimport)
#endif

#define LIB_drawNBNfun_CPP_API PUBLIC_drawNBNfun_CPP_API

#else

#if !defined(LIB_drawNBNfun_CPP_API)
#if defined(LIB_drawNBNfun_C_API)
#define LIB_drawNBNfun_CPP_API LIB_drawNBNfun_C_API
#else
#define LIB_drawNBNfun_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_drawNBNfun_CPP_API void MW_CALL_CONV drawNBNfun(const mwArray& figureId, const mwArray& edge_a, const mwArray& edge_b, const mwArray& x, const mwArray& y, const mwArray& z, const mwArray& node_fit, const mwArray& numColor, const mwArray& basicNodeSize);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
