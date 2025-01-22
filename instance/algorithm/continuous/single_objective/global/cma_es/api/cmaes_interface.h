/* --------------------------------------------------------- */
/* --- File: cmaes_interface.h - Author: Nikolaus Hansen --- */
/* ---------------------- last modified:  IV 2007        --- */
/* --------------------------------- by: Nikolaus Hansen --- */
/* --------------------------------------------------------- */
/*   
     CMA-ES for non-linear function minimization. 

     Copyright (C) 1996, 2003, 2007 Nikolaus Hansen. 
     e-mail: hansen AT lri.fr
     
     Documentation: see file docfunctions.txt
     
     License: see file cmaes.c
*/
#include "cmaes.h"
#include "../../../../../../../core/random/newran.h"

/* --------------------------------------------------------- */
/* ------------------ Interface ---------------------------- */
/* --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

/* --- initialization, constructors, destructors --- */
double* cmaes_init(cmaes_t *, int dimension , double *xstart, 
		double *stddev/*, long seed*/, int lambda,
		const char *input_parameter_filename, ofec::Random *rnd);
void cmaes_init_para(cmaes_t *, int dimension , double *xstart, 
		double *stddev/*, long seed*/, int lambda,
		const char *input_parameter_filename);
double* cmaes_init_final(cmaes_t *, ofec::Random *rnd);
void cmaes_resumeDistribution(cmaes_t *evo_ptr, char *filename);
void cmaes_exit(cmaes_t *);

/* --- core functions --- */
double * const * cmaes_SamplePopulation(cmaes_t *, ofec::Random *rnd);
double *         cmaes_UpdateDistribution(cmaes_t *, 
					  const double *rgFitnessValues);
const char *     cmaes_TestForTermination(cmaes_t *);

/* --- additional functions --- */
double * const * cmaes_ReSampleSingle( cmaes_t *t, int index, ofec::Random *rnd);
double const *   cmaes_ReSampleSingle_old(cmaes_t *, double *rgx, ofec::Random *rnd);
double *         cmaes_SampleSingleInto( cmaes_t *t, double *rgx, ofec::Random *rnd);
void             cmaes_UpdateEigensystem(cmaes_t *, int flgforce);
void             cmaes_ResizeLambda(cmaes_t* t, int lambda);
void             cmaes_InitDistribution(cmaes_t* t, const double* rgFunVal);
bool             cmaes_CheckEqualFitness(cmaes_t *t, const double *rgFunVal);

/* --- getter functions --- */
double         cmaes_Get(cmaes_t *, char const *keyword);
const double * cmaes_GetPtr(cmaes_t *, char const *keyword); /* e.g. "xbestever" */
double *       cmaes_GetNew( cmaes_t *t, char const *keyword); /* user is responsible to free */
double *       cmaes_GetInto( cmaes_t *t, char const *keyword, double *mem); /* allocs if mem==NULL, user is responsible to free */
double* const* cmaes_GetPop(cmaes_t* t);
double*        cmaes_GetArFunvals(cmaes_t* t);

/* --- online control and output --- */
void           cmaes_ReadSignals(cmaes_t *, char const *filename);
void           cmaes_WriteToFile(cmaes_t *, const char *szKeyWord,
                                 const char *output_filename); 
char *         cmaes_SayHello(cmaes_t *);
/* --- misc --- */
double *       cmaes_NewDouble(int n); /* user is responsible to free */
void           cmaes_FATAL(char const *s1, char const *s2, char const *s3, 
			   char const *s4);

#ifdef __cplusplus
} // end extern "C"
#endif

