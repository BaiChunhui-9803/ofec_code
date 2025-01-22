//
// MATLAB Compiler: 8.2 (R2021a)
// Date: Sat Jun  3 14:30:06 2023
// Arguments:
// "-B""macro_default""-W""cpplib:drawdf,all,version=1.0""-T""link:lib""-d""E:\D
// iao_Yiya\code\OFEC\utility\matlab\instance\drawDF\drawdf\for_testing""-v""\\1
// 72.24.242.8\share\Student\2018\YiyaDiao\nbn_code\matlab_visual\drawdf.m"
//

#define EXPORTING_drawdf 1
#include "drawdf.h"

static HMCRINSTANCE _mcr_inst = NULL; /* don't use nullptr; this may be either C or C++ */

#if defined( _MSC_VER) || defined(__LCC__) || defined(__MINGW64__)
#ifdef __LCC__
#undef EXTERN_C
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#define NOMINMAX
#include <windows.h>
#undef interface

static char path_to_dll[_MAX_PATH];

BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, void *pv)
{
    if (dwReason == DLL_PROCESS_ATTACH)
    {
        if (GetModuleFileName(hInstance, path_to_dll, _MAX_PATH) == 0)
            return FALSE;
    }
    else if (dwReason == DLL_PROCESS_DETACH)
    {
    }
    return TRUE;
}
#endif
#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

static int mclDefaultPrintHandler(const char *s)
{
    return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern C block */
#endif

#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

static int mclDefaultErrorHandler(const char *s)
{
    int written = 0;
    size_t len = 0;
    len = strlen(s);
    written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
    if (len > 0 && s[ len-1 ] != '\n')
        written += mclWrite(2 /* stderr */, "\n", sizeof(char));
    return written;
}

#ifdef __cplusplus
} /* End extern C block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_drawdf_C_API
#define LIB_drawdf_C_API /* No special import/export declaration */
#endif

LIB_drawdf_C_API 
bool MW_CALL_CONV drawdfInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler)
{
    int bResult = 0;
    if (_mcr_inst)
        return true;
    if (!mclmcrInitialize())
        return false;
    if (!GetModuleFileName(GetModuleHandle("drawdf"), path_to_dll, _MAX_PATH))
        return false;
    {
        mclCtfStream ctfStream = 
            mclGetEmbeddedCtfStream(path_to_dll);
        if (ctfStream) {
            bResult = mclInitializeComponentInstanceEmbedded(&_mcr_inst,
                                                             error_handler, 
                                                             print_handler,
                                                             ctfStream);
            mclDestroyStream(ctfStream);
        } else {
            bResult = 0;
        }
    }  
    if (!bResult)
    return false;
    return true;
}

LIB_drawdf_C_API 
bool MW_CALL_CONV drawdfInitialize(void)
{
    return drawdfInitializeWithHandlers(mclDefaultErrorHandler, mclDefaultPrintHandler);
}

LIB_drawdf_C_API 
void MW_CALL_CONV drawdfTerminate(void)
{
    if (_mcr_inst)
        mclTerminateInstance(&_mcr_inst);
}

LIB_drawdf_C_API 
void MW_CALL_CONV drawdfPrintStackTrace(void) 
{
    char** stackTrace;
    int stackDepth = mclGetStackTrace(&stackTrace);
    int i;
    for(i=0; i<stackDepth; i++)
    {
        mclWrite(2 /* stderr */, stackTrace[i], sizeof(char)*strlen(stackTrace[i]));
        mclWrite(2 /* stderr */, "\n", sizeof(char)*strlen("\n"));
    }
    mclFreeStackTrace(&stackTrace, stackDepth);
}


LIB_drawdf_C_API 
bool MW_CALL_CONV mlxDrawdf(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    return mclFeval(_mcr_inst, "drawdf", nlhs, plhs, nrhs, prhs);
}

LIB_drawdf_CPP_API 
void MW_CALL_CONV drawdf(const mwArray& figureId, const mwArray& df_mat, const mwArray& 
                         numColor, const mwArray& basicNodeSize)
{
    mclcppMlfFeval(_mcr_inst, "drawdf", 0, 0, 4, &figureId, &df_mat, &numColor, &basicNodeSize);
}

