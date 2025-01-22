//
// MATLAB Compiler: 8.2 (R2021a)
// Date: Mon Oct  9 20:22:08 2023
// Arguments:
// "-B""macro_default""-W""cpplib:drawSaveNBNInMesh2Imagefun,all,version=1.0""-T
// ""link:lib""-d""E:\Diao_Yiya\code\OFEC\utility\matlab\instance\showNBN\drawSa
// veNBNMesh\drawSaveNBNInMesh2Imagefun\for_testing""-v""E:\Diao_Yiya\code\OFEC\
// utility\matlab\script\nbn_figure_out\drawSaveNBNInMesh2Imagefun.m"
//

#ifndef drawSaveNBNInMesh2Imagefun_h
#define drawSaveNBNInMesh2Imagefun_h 1

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
#ifndef LIB_drawSaveNBNInMesh2Imagefun_C_API 
#define LIB_drawSaveNBNInMesh2Imagefun_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_drawSaveNBNInMesh2Imagefun_C_API 
bool MW_CALL_CONV drawSaveNBNInMesh2ImagefunInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_drawSaveNBNInMesh2Imagefun_C_API 
bool MW_CALL_CONV drawSaveNBNInMesh2ImagefunInitialize(void);

extern LIB_drawSaveNBNInMesh2Imagefun_C_API 
void MW_CALL_CONV drawSaveNBNInMesh2ImagefunTerminate(void);

extern LIB_drawSaveNBNInMesh2Imagefun_C_API 
void MW_CALL_CONV drawSaveNBNInMesh2ImagefunPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_drawSaveNBNInMesh2Imagefun_C_API 
bool MW_CALL_CONV mlxDrawSaveNBNInMesh2Imagefun(int nlhs, mxArray *plhs[], int nrhs, 
                                                mxArray *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_drawSaveNBNInMesh2Imagefun
#define PUBLIC_drawSaveNBNInMesh2Imagefun_CPP_API __declspec(dllexport)
#else
#define PUBLIC_drawSaveNBNInMesh2Imagefun_CPP_API __declspec(dllimport)
#endif

#define LIB_drawSaveNBNInMesh2Imagefun_CPP_API PUBLIC_drawSaveNBNInMesh2Imagefun_CPP_API

#else

#if !defined(LIB_drawSaveNBNInMesh2Imagefun_CPP_API)
#if defined(LIB_drawSaveNBNInMesh2Imagefun_C_API)
#define LIB_drawSaveNBNInMesh2Imagefun_CPP_API LIB_drawSaveNBNInMesh2Imagefun_C_API
#else
#define LIB_drawSaveNBNInMesh2Imagefun_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_drawSaveNBNInMesh2Imagefun_CPP_API void MW_CALL_CONV drawSaveNBNInMesh2Imagefun(const mwArray& figureId, const mwArray& filepath, const mwArray& filetype, const mwArray& edge_a, const mwArray& edge_b, const mwArray& x, const mwArray& y, const mwArray& z, const mwArray& node_fit, const mwArray& numColor, const mwArray& basicNodeSize, const mwArray& linewidth, const mwArray& mesh_x, const mwArray& mesh_y, const mwArray& mesh_z, const mwArray& conT);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
