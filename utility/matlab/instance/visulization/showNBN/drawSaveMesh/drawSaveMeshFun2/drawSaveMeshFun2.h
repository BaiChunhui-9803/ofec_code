//
// MATLAB Compiler: 8.2 (R2021a)
// Date: Mon Aug 28 23:20:03 2023
// Arguments:
// "-B""macro_default""-W""cpplib:drawSaveMeshFun2,all,version=1.0""-T""link:lib
// ""-d""E:\Diao_Yiya\code\OFEC\utility\matlab\instance\showNBN\drawSaveMesh\dra
// wSaveMeshFun2\for_testing""-v""E:\Diao_Yiya\code\OFEC\utility\matlab\script\n
// bn_figure_out\drawSaveMeshFun2.m"
//

#ifndef drawSaveMeshFun2_h
#define drawSaveMeshFun2_h 1

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
#ifndef LIB_drawSaveMeshFun2_C_API 
#define LIB_drawSaveMeshFun2_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_drawSaveMeshFun2_C_API 
bool MW_CALL_CONV drawSaveMeshFun2InitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_drawSaveMeshFun2_C_API 
bool MW_CALL_CONV drawSaveMeshFun2Initialize(void);

extern LIB_drawSaveMeshFun2_C_API 
void MW_CALL_CONV drawSaveMeshFun2Terminate(void);

extern LIB_drawSaveMeshFun2_C_API 
void MW_CALL_CONV drawSaveMeshFun2PrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_drawSaveMeshFun2_C_API 
bool MW_CALL_CONV mlxDrawSaveMeshFun2(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                      *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_drawSaveMeshFun2
#define PUBLIC_drawSaveMeshFun2_CPP_API __declspec(dllexport)
#else
#define PUBLIC_drawSaveMeshFun2_CPP_API __declspec(dllimport)
#endif

#define LIB_drawSaveMeshFun2_CPP_API PUBLIC_drawSaveMeshFun2_CPP_API

#else

#if !defined(LIB_drawSaveMeshFun2_CPP_API)
#if defined(LIB_drawSaveMeshFun2_C_API)
#define LIB_drawSaveMeshFun2_CPP_API LIB_drawSaveMeshFun2_C_API
#else
#define LIB_drawSaveMeshFun2_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_drawSaveMeshFun2_CPP_API void MW_CALL_CONV drawSaveMeshFun2(const mwArray& x, const mwArray& y, const mwArray& z, const mwArray& conT, const mwArray& figureId, const mwArray& filepath);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
