//
// MATLAB Compiler: 8.2 (R2021a)
// Date: Mon Oct  9 18:53:52 2023
// Arguments:
// "-B""macro_default""-W""cpplib:drawSaveNBN2Imagefun3,all,version=1.0""-T""lin
// k:lib""-d""E:\Diao_Yiya\code\OFEC\utility\matlab\instance\showNBN\drawSaveNBN
// 3\drawSaveNBN2Imagefun3\for_testing""-v""E:\Diao_Yiya\code\OFEC\utility\matla
// b\script\nbn_figure_out\drawSaveNBN2Imagefun3.m"
//

#define EXPORTING_drawSaveNBN2Imagefun3 1
#include "drawSaveNBN2Imagefun3.h"

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
#ifndef LIB_drawSaveNBN2Imagefun3_C_API
#define LIB_drawSaveNBN2Imagefun3_C_API /* No special import/export declaration */
#endif

LIB_drawSaveNBN2Imagefun3_C_API 
bool MW_CALL_CONV drawSaveNBN2Imagefun3InitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler)
{
    int bResult = 0;
    if (_mcr_inst)
        return true;
    if (!mclmcrInitialize())
        return false;
    if (!GetModuleFileName(GetModuleHandle("drawSaveNBN2Imagefun3"), path_to_dll, _MAX_PATH))
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

LIB_drawSaveNBN2Imagefun3_C_API 
bool MW_CALL_CONV drawSaveNBN2Imagefun3Initialize(void)
{
    return drawSaveNBN2Imagefun3InitializeWithHandlers(mclDefaultErrorHandler, 
                                                     mclDefaultPrintHandler);
}

LIB_drawSaveNBN2Imagefun3_C_API 
void MW_CALL_CONV drawSaveNBN2Imagefun3Terminate(void)
{
    if (_mcr_inst)
        mclTerminateInstance(&_mcr_inst);
}

LIB_drawSaveNBN2Imagefun3_C_API 
void MW_CALL_CONV drawSaveNBN2Imagefun3PrintStackTrace(void) 
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


LIB_drawSaveNBN2Imagefun3_C_API 
bool MW_CALL_CONV mlxDrawSaveNBN2Imagefun3(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                           *prhs[])
{
    return mclFeval(_mcr_inst, "drawSaveNBN2Imagefun3", nlhs, plhs, nrhs, prhs);
}

LIB_drawSaveNBN2Imagefun3_CPP_API 
void MW_CALL_CONV drawSaveNBN2Imagefun3(const mwArray& figureId, const mwArray& filepath, 
                                        const mwArray& filetype, const mwArray& edge_a, 
                                        const mwArray& edge_b, const mwArray& x, const 
                                        mwArray& y, const mwArray& z, const mwArray& 
                                        node_fit, const mwArray& numColor, const mwArray& 
                                        basicNodeSize, const mwArray& linewidth)
{
    mclcppMlfFeval(_mcr_inst, "drawSaveNBN2Imagefun3", 0, 0, 12, &figureId, &filepath, &filetype, &edge_a, &edge_b, &x, &y, &z, &node_fit, &numColor, &basicNodeSize, &linewidth);
}

