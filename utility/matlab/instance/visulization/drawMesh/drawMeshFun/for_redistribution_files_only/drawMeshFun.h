//
// MATLAB Compiler: 8.2 (R2021a)
// Date: Mon May 29 20:48:33 2023
// Arguments:
// "-B""macro_default""-W""cpplib:drawMeshFun,all,version=1.0""-T""link:lib""-d"
// "E:\Diao_Yiya\code\OFEC\utility\matlab\instance\drawMesh\drawMeshFun\for_test
// ing""-v""\\172.24.242.8\share\Student\2018\YiyaDiao\nbn_code\matlab_visual\dr
// awMeshFun.m"
//

#ifndef drawMeshFun_h
#define drawMeshFun_h 1

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
#ifndef LIB_drawMeshFun_C_API 
#define LIB_drawMeshFun_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_drawMeshFun_C_API 
bool MW_CALL_CONV drawMeshFunInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_drawMeshFun_C_API 
bool MW_CALL_CONV drawMeshFunInitialize(void);

extern LIB_drawMeshFun_C_API 
void MW_CALL_CONV drawMeshFunTerminate(void);

extern LIB_drawMeshFun_C_API 
void MW_CALL_CONV drawMeshFunPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_drawMeshFun_C_API 
bool MW_CALL_CONV mlxDrawMeshFun(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_drawMeshFun
#define PUBLIC_drawMeshFun_CPP_API __declspec(dllexport)
#else
#define PUBLIC_drawMeshFun_CPP_API __declspec(dllimport)
#endif

#define LIB_drawMeshFun_CPP_API PUBLIC_drawMeshFun_CPP_API

#else

#if !defined(LIB_drawMeshFun_CPP_API)
#if defined(LIB_drawMeshFun_C_API)
#define LIB_drawMeshFun_CPP_API LIB_drawMeshFun_C_API
#else
#define LIB_drawMeshFun_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_drawMeshFun_CPP_API void MW_CALL_CONV drawMeshFun(const mwArray& x, const mwArray& y, const mwArray& z, const mwArray& conT);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
