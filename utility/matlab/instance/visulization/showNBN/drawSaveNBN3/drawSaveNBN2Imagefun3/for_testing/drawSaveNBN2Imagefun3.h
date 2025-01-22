//
// MATLAB Compiler: 8.2 (R2021a)
// Date: Mon Oct  9 18:53:52 2023
// Arguments:
// "-B""macro_default""-W""cpplib:drawSaveNBN2Imagefun3,all,version=1.0""-T""lin
// k:lib""-d""E:\Diao_Yiya\code\OFEC\utility\matlab\instance\showNBN\drawSaveNBN
// 3\drawSaveNBN2Imagefun3\for_testing""-v""E:\Diao_Yiya\code\OFEC\utility\matla
// b\script\nbn_figure_out\drawSaveNBN2Imagefun3.m"
//

#ifndef drawSaveNBN2Imagefun3_h
#define drawSaveNBN2Imagefun3_h 1

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
#ifndef LIB_drawSaveNBN2Imagefun3_C_API 
#define LIB_drawSaveNBN2Imagefun3_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_drawSaveNBN2Imagefun3_C_API 
bool MW_CALL_CONV drawSaveNBN2Imagefun3InitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_drawSaveNBN2Imagefun3_C_API 
bool MW_CALL_CONV drawSaveNBN2Imagefun3Initialize(void);

extern LIB_drawSaveNBN2Imagefun3_C_API 
void MW_CALL_CONV drawSaveNBN2Imagefun3Terminate(void);

extern LIB_drawSaveNBN2Imagefun3_C_API 
void MW_CALL_CONV drawSaveNBN2Imagefun3PrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_drawSaveNBN2Imagefun3_C_API 
bool MW_CALL_CONV mlxDrawSaveNBN2Imagefun3(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                           *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_drawSaveNBN2Imagefun3
#define PUBLIC_drawSaveNBN2Imagefun3_CPP_API __declspec(dllexport)
#else
#define PUBLIC_drawSaveNBN2Imagefun3_CPP_API __declspec(dllimport)
#endif

#define LIB_drawSaveNBN2Imagefun3_CPP_API PUBLIC_drawSaveNBN2Imagefun3_CPP_API

#else

#if !defined(LIB_drawSaveNBN2Imagefun3_CPP_API)
#if defined(LIB_drawSaveNBN2Imagefun3_C_API)
#define LIB_drawSaveNBN2Imagefun3_CPP_API LIB_drawSaveNBN2Imagefun3_C_API
#else
#define LIB_drawSaveNBN2Imagefun3_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_drawSaveNBN2Imagefun3_CPP_API void MW_CALL_CONV drawSaveNBN2Imagefun3(const mwArray& figureId, const mwArray& filepath, const mwArray& filetype, const mwArray& edge_a, const mwArray& edge_b, const mwArray& x, const mwArray& y, const mwArray& z, const mwArray& node_fit, const mwArray& numColor, const mwArray& basicNodeSize, const mwArray& linewidth);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
