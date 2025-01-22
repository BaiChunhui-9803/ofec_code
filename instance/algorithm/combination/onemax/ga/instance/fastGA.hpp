/********* Begin Register Information **********
{
    "name": "FastGA-OneMax",
    "identifier": "FastGA_OneMax",
    "problem tags": [ "ComOP", "GOP", "SOP", "MMOP", "OneMax" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*-------------------------------------------------------------------------------
*
*
*-------------------------------------------------------------------------------
*
*
*************************************************************************/


#include "../ga_onemax.h"
#include "../../../../../../utility/typevar/typevar.h"
namespace ofec
{

  
    class FastGA_OneMax : public GeneticAlgorithmOneMax
    {

    public:

        static void setParameter(ParamMap& param, int mu, int lambda)
      {

        param["mu"] = mu;
        param["lambda"] = lambda;

      }


        virtual void updateParamTypes(paramTypeMap& param) {
            // VisualParamterBase* curPara;
            {
                VisualParamterIntSlider curPara;
                curPara.m_key = "mu";
                curPara.setRange(1, 30, 1, 1);
                param.emplace_back(new VisualParamterIntSlider(curPara));
            }

            {
                VisualParamterIntSlider curPara;
                curPara.m_key = "lambda";
                curPara.setRange(1, 30, 1, 1);
                param.emplace_back(new VisualParamterIntSlider(curPara));
            }

            //{
            //    VisualParamterRealSlider curPara;
            //    curPara.m_key = "mutation_rate_scale_";
            //    curPara.setRange(1, 30, 1, 1);
            //    param.emplace_back(new VisualParamterRealSlider(curPara));
            //}
        }

        void initialize_() override {
            GeneticAlgorithmOneMax::initialize_();
            int mu;
            int lambda;
            //double mutation_rate_scale;
            auto& v = *m_param;
            mu = v.get<int>("mu");
            lambda = v.get<int>("lambda");
         //   mutation_rate_scale = v.get<double>("mutation_rate_scale_");
            this->set_mu(mu);
            this->set_lambda(lambda);
            this->set_crossover_mutation_r("OR");
            this->set_crossover_probability(0);
            this->set_mutation_operator("POWERLAWSAMPLE");
            this->set_beta_f(1.5);
            this->set_selection_operator("BESTPLUS");
        }

    };
  
}
