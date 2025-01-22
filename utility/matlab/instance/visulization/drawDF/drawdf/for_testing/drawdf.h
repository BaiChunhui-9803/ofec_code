//
// MATLAB Compiler: 8.2 (R2021a)
// Date: Sat Jun  3 14:30:06 2023
// Arguments:
// "-B""macro_default""-W""cpplib:drawdf,all,version=1.0""-T""link:lib""-d""E:\D
// iao_Yiya\code\OFEC\utility\matlab\instance\drawDF\drawdf\for_testing""-v""\\1
// 72.24.242.8\share\Student\2018\YiyaDiao\nbn_code\matlab_visual\drawdf.m"
//

#ifndef drawdf_h
#define drawdf_h 1

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
#ifndef LIB_drawdf_C_API 
#define LIB_drawdf_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_drawdf_C_API 
bool MW_CALL_CONV drawdfInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_drawdf_C_API 
bool MW_CALL_CONV drawdfInitialize(void);

extern LIB_drawdf_C_API 
void MW_CALL_CONV drawdfTerminate(void);

extern LIB_drawdf_C_API 
void MW_CALL_CONV drawdfPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_drawdf_C_API 
bool MW_CALL_CONV mlxDrawdf(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_drawdf
#define PUBLIC_drawdf_CPP_API __declspec(dllexport)
#else
#define PUBLIC_drawdf_CPP_API __declspec(dllimport)
#endif

#define LIB_drawdf_CPP_API PUBLIC_drawdf_CPP_API

#else

#if !defined(LIB_drawdf_CPP_API)
#if defined(LIB_drawdf_C_API)
#define LIB_drawdf_CPP_API LIB_drawdf_C_API
#else
#define LIB_drawdf_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_drawdf_CPP_API void MW_CALL_CONV drawdf(const mwArray& figureId, const mwArray& df_mat, const mwArray& numColor, const mwArray& basicNodeSize);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
