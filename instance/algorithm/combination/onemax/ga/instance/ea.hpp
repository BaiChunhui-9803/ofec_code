/********* Begin Register Information **********
{
    "name": "staticEA-OneMax",
    "identifier": "staticEA_OneMax",
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

#pragma once
//#include "../geneticAlgorithm.hpp"
#include "../ga_onemax.h"
#include "../../../../../../utility/typevar/typevar.h"
namespace ofec {


    class staticEA_OneMax : public GeneticAlgorithmOneMax
    {

    public:
        static void setParameter(ParamMap& param, int mu, int lambda, double mutation_rate_scale = 1)
        {
            param["mu"] = mu;
            param["lambda"] = lambda;
            param["mutation_rate_scale_"] = mutation_rate_scale;
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

            {
                VisualParamterRealSlider curPara;
                curPara.m_key = "mutation_rate_scale_";
                curPara.setRange(1, 30, 1, 0.1);
                param.emplace_back(new VisualParamterRealSlider(curPara));
            }
        }

        void initialize_() override {
            GeneticAlgorithmOneMax::initialize_();



            int mu;
            int lambda;
            double mutation_rate_scale;
            auto& v = *m_param;
            mu = v.get<int>("mu");
            lambda = v.get<int>("lambda");
            mutation_rate_scale = v.get<double>("mutation_rate_scale_");
            GeneticAlgorithmOneMax::set_mu(mu);
            this->set_lambda(lambda);
            this->set_crossover_mutation_r("OR");
            this->set_crossover_probability(0);
            this->set_mutation_operator("BINOMIALSAMPLE");
            this->set_selection_operator("BESTPLUS");
            //mutation_rate_scale_ = mutation_rate_scale;

            this->set_mutation_rate(mutation_rate_scale / static_cast<double>(m_problem->numVariables()));

        }


    };
  
}
