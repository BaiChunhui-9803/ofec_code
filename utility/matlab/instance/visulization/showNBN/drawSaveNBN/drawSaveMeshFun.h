//
// MATLAB Compiler: 8.2 (R2021a)
// Date: Mon Aug 28 21:54:02 2023
// Arguments:
// "-B""macro_default""-W""cpplib:drawSaveMeshFun,all,version=1.0""-T""link:lib"
// "-d""E:\Diao_Yiya\code\OFEC\utility\matlab\instance\showNBN\drawSaveNBN\drawS
// aveMeshFun\for_testing""-v""E:\Diao_Yiya\code\OFEC\utility\matlab\script\stn_
// visualization\drawSaveMeshFun.m"
//

#ifndef drawSaveMeshFun_h
#define drawSaveMeshFun_h 1

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
#ifndef LIB_drawSaveMeshFun_C_API 
#define LIB_drawSaveMeshFun_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_drawSaveMeshFun_C_API 
bool MW_CALL_CONV drawSaveMeshFunInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_drawSaveMeshFun_C_API 
bool MW_CALL_CONV drawSaveMeshFunInitialize(void);

extern LIB_drawSaveMeshFun_C_API 
void MW_CALL_CONV drawSaveMeshFunTerminate(void);

extern LIB_drawSaveMeshFun_C_API 
void MW_CALL_CONV drawSaveMeshFunPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_drawSaveMeshFun_C_API 
bool MW_CALL_CONV mlxDrawSaveMeshFun(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                     *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_drawSaveMeshFun
#define PUBLIC_drawSaveMeshFun_CPP_API __declspec(dllexport)
#else
#define PUBLIC_drawSaveMeshFun_CPP_API __declspec(dllimport)
#endif

#define LIB_drawSaveMeshFun_CPP_API PUBLIC_drawSaveMeshFun_CPP_API

#else

#if !defined(LIB_drawSaveMeshFun_CPP_API)
#if defined(LIB_drawSaveMeshFun_C_API)
#define LIB_drawSaveMeshFun_CPP_API LIB_drawSaveMeshFun_C_API
#else
#define LIB_drawSaveMeshFun_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_drawSaveMeshFun_CPP_API void MW_CALL_CONV drawSaveMeshFun(const mwArray& meshData, const mwArray& conT, const mwArray& figureId, const mwArray& filepath);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
