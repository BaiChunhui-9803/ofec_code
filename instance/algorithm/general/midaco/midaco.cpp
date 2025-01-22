// Include directory "api"

#include "midaco.h"
#include "../../../../core/problem/continuous/continuous.h"


#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

extern"C" { int midaco(long int *, long int *, long int *, long int *, long int *,
    long int *, double *, double *, double *, double *, double *,
    long int *, long int *, double *, double *, long int *,
    long int *, long int *, double *, long int *, char *); }
/***********************************************************************/
extern"C" { int midaco_print(int, long int, long int, long int *, long int *, double *,
    double *, double *, double *, double *, long int, long int,
    long int, long int, long int, double *, double *,
    long int, long int, double *, long int, char *); }

namespace ofec {
	void MIDACO::initialize_(Environment *env) {
		Algorithm::initialize_(env);
	}

	void MIDACO::run_(Environment *env) {
        m_sol.reset(m_problem->createSolution());
        m_sol->initialize(env, m_random.get());

        /* Variable and Workspace Declarations */
        long int o, n, ni, m, me, maxeval, maxtime, printeval, save2file, iflag, istop;
        long int liw, lrw, lpf, i, iw[5000], p = 1; double rw[20000], pf[20000];
        double   f[10], g[1000], x[1000], xl[1000], xu[1000], param[13];
        char key[] = "MIDACO_LIMITED_VERSION___[CREATIVE_COMMONS_BY-NC-ND_LICENSE]";

        /* STEP 1.A: Problem dimensions
        ******************************/
        o = m_problem->numberObjectives(); /* Number of objectives                          */
        n = env->problem()->numberVariables(); /* Number of variables (in total)                */
        ni = 0; /* Number of integer variables (0 <= ni <= n)    */
        m = m_problem->numberConstraints(); /* Number of constraints (in total)              */
        me = 0; /* Number of equality constraints (0 <= me <= m) */
        for (i = 0; i < m; ++i) {
            if (m_problem->isEqualityConstraint(i)) {
                me++;
            }
        }

        /* STEP 1.B: Lower and upper bounds 'xl' & 'xu'
        **********************************************/
        for (i = 0; i < n; i++) {
            xl[i] = CAST_CONOP(env->problem())->domain().range(i).limit.first;
            xu[i] = CAST_CONOP(env->problem())->domain().range(i).limit.second;
        }

        /* STEP 1.C: Starting point 'x'
        ******************************/
        for (i = 0; i < n; i++) {
            x[i] = dynamic_cast<VariableVector<>&>(m_sol->variableBase())[i]; /* Here for example: starting point = lower bounds */
        }

        /* STEP 2.A: Stopping criteria
        *****************************/
        maxeval = m_maximum_evalutions;    /* Maximum number of function evaluation (e.g. 1000000)  */
        maxtime = 60 * 60 * 24; /* Maximum time limit in Seconds (e.g. 1 Day = 60*60*24) */

        /* STEP 2.B: Printing options
****************************/
        printeval = 1000; /* Print-Frequency for current best solution (e.g. 1000) */
        save2file = 0;    /* Save SCREEN and SOLUTION to TXT-files [ 0=NO/ 1=YES]  */

        /*****************************************************************/
        /***  Step 3: Choose MIDACO parameters (FOR ADVANCED USERS)    ***/
        /*****************************************************************/

        param[0] = 0.0;  /* ACCURACY  */
        param[1] = m_random->seed;  /* SEED      */
        param[2] = 0.0;  /* FSTOP     */
        param[3] = 0.0;  /* ALGOSTOP  */
        param[4] = 0.0;  /* EVALSTOP  */
        param[5] = 0.0;  /* FOCUS     */
        param[6] = 0.0;  /* ANTS      */
        param[7] = 0.0;  /* KERNEL    */
        param[8] = 0.0;  /* ORACLE    */
        param[9] = 0.0;  /* PARETOMAX */
        param[10] = 0.0;  /* EPSILON   */
        param[11] = 0.0;  /* BALANCE   */
        param[12] = 0.0;  /* CHARACTER */

        std::vector<Real> pos(o);
        for (size_t i = 0; i < o; ++i) {
            pos[i] = m_problem->optimizeMode(i) == OptimizeMode::kMinimize ? 1 : -1;
        }

        /* Workspace length calculation */
        lrw = sizeof(rw) / sizeof(double);
        lpf = sizeof(pf) / sizeof(double);
        liw = sizeof(iw) / sizeof(long int);

        /* Print midaco headline and basic information */
        midaco_print(1, printeval, save2file, &iflag, &istop, &*f, &*g, &*x, &*xl, &*xu,
            o, n, ni, m, me, &*rw, &*pf, maxeval, maxtime, &*param, p, &*key);

#ifdef OFEC_DATUM_MULTI_POP_H
        datum::g_multi_pop.clear();
        datum::g_multi_pop.resize(1);
        datum::g_multi_pop[0].push_back(m_sol.get());
#endif

		while (!terminating() && istop == 0) { /*~~~ Start of the reverse communication loop ~~~*/
            /* Evaluate objective function */
			m_sol->evaluate(env);
            datumUpdated(env);
            for (size_t i = 0; i < o; ++i) {
                f[i] = pos[i] * m_sol->objective(i);
            }
            for (size_t i = 0; i < m; ++i) {
                g[i] = m_sol->constraint(i);
            }
            /* Call MIDACO */
            midaco(&p, &o, &n, &ni, &m, &me, &*x, &*f, &*g, &*xl, &*xu, &iflag,
                &istop, &*param, &*rw, &lrw, &*iw, &liw, &*pf, &lpf, &*key);
            for (size_t i = 0; i < n; ++i) {
                dynamic_cast<VariableVector<>&>(m_sol->variableBase())[i] = x[i];
            }
            /* Call MIDACO printing routine */
            midaco_print(2, printeval, save2file, &iflag, &istop, &*f, &*g, &*x, &*xl, &*xu,
                o, n, ni, m, me, &*rw, &*pf, maxeval, maxtime, &*param, p, &*key);
		} /*~~~End of the reverse communication loop ~~~*/ 
	}
}