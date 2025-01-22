#include "nl_shade_rsp_mid.h"
#include "optimizer.h"
#include "../../../../../../core/problem/continuous/continuous.h"


namespace ofec {


    void NL_SHADE_RSP_MID::initialize_(Environment* env)  {
        Algorithm::initialize_(env);
        ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
        auto dimensionality = con_pro->numberVariables();

        POP_SIZE = dimensionality * 5;
        oldPopSize = POP_SIZE;
    }
    void NL_SHADE_RSP_MID::run_(Environment* env) {

        ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
        auto dimensionality = con_pro->numberVariables();

        Optimizer OptZ;
        ////////////////popSize     NVars  func     Run  bench     memory    arch size

#if defined (SR_SAVE_AS_NEAREST_RESTART) || defined (K_MEANS_AS_NEAREST)

        OptZ.Initialize(POP_SIZE, dimensionality, 20 * dimensionality, 2.1, m_random.get());
        //   isRestart = false;
        while (!terminating()) {
            OptZ.MainCycle(fopt, env, m_random.get());

            POP_SIZE = 400; //best

            OptZ.restart(POP_SIZE, dimensionality, 20 * dimensionality, 2.1, env, m_random.get());
            //  isRestart = true;

        };
        POP_SIZE = oldPopSize;
        //    isRestart = false;
#endif                    

        OptZ.Clean();
    }
}