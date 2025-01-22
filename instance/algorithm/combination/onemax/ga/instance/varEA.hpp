/********* Begin Register Information **********
{
    "name": "VarEA-OneMax",
    "identifier": "VarEA_OneMax",
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
#include "../ga_onemax.h"
namespace ofec {
    //namespace modularGA
    //{

    //    namespace ga
    //    {
            class VarEA_OneMax : public GeneticAlgorithmOneMax
            {
                double init_F = 0.98;
                std::vector<int> best_ind, x;
                double r_n = this->get_r_n();
                double best_f_x, f_x;
                int flip_n, best_r_n;
                double F = init_F;
                int c = 0;
            public:


                static void setParameter(ParamMap& param, 
                    int lambda, double init_r_n_ = 1.0)
                {
                    param["lambda"] = lambda;
                    param["init_r_n_"] = init_r_n_;
                }
                virtual void updateParamTypes(paramTypeMap& param) {
                    // VisualParamterBase* curPara;
                    //{
                    //    VisualParamterIntSlider curPara;
                    //    curPara.m_key = "mu";
                    //    curPara.setRange(1, 30, 1, 1);
                    //    param.emplace_back(new VisualParamterIntSlider(curPara));
                    //}

                    {
                        VisualParamterIntSlider curPara;
                        curPara.m_key = "lambda";
                        curPara.setRange(1, 30, 1, 1);
                        param.emplace_back(new VisualParamterIntSlider(curPara));
                    }

                    {
                        VisualParamterRealSlider curPara;
                        curPara.m_key = "init_r_n_";
                        curPara.setRange(1, 30, 1, 0.1);
                        param.emplace_back(new VisualParamterRealSlider(curPara));
                    }
                }

                void initialize_() override {
                    GeneticAlgorithmOneMax::initialize_();


                    int lambda;
                    double init_r_n_;
                    auto& v = *m_param;
                    lambda = v.get<int>("lambda");
                    init_r_n_ = v.get<double>("init_r_n_");


                    this->set_mu(1);
                    this->set_lambda(lambda);
                    this->set_mutation_operator("NORMALSAMPLE");
                    this->set_r_n(init_r_n_);

                    init_F = 0.98;
                    //std::vector<int> best_ind, x;
                    r_n = this->get_r_n();
 /*                   double best_f_x, f_x;
                    int flip_n, best_r_n;*/
                    F = init_F;
                    c = 0;

                    this->Preparation(m_problem->numVariables());
                    this->Initialization();




                }

                virtual void evolve() override
                {
 
 


                   // while (!this->Termination())
                    {
                        this->update_generation();
                        best_r_n = r_n;
                        best_f_x = std::numeric_limits<int>::min();
                        this->set_r_n(r_n);
                        this->set_sigma_n(sqrt(pow(F, c) *
                            static_cast<double>(r_n) *
                            (1.0 - r_n / static_cast<double>(m_problem->numVariables()))));
                        for (size_t i = 0; i < this->get_lambda(); ++i)
                        {
                            x = this->get_parents_population()[0];
                            flip_n = this->DoMutation(x,m_random.get());
                            f_x = this->Evaluate(x);
                            if (f_x > best_f_x)
                            {
                                best_f_x = f_x;
                                r_n = flip_n;
                                best_ind = x;
                            }

               /*             if (this->Termination())
                                break;*/
                        }

      /*                  if (this->Termination())
                            break;*/

                        if (best_r_n == r_n)
                        {
                            c++;
                        }
                        else
                        {
                            c = 0;
                        }

                        if (best_f_x >= this->get_parents_fitness()[0])
                        {
                            this->set_parents_population(best_ind, 0);
                            this->set_parents_fitness(best_f_x, 0);
                        }
                    }
                }

            };
    //    }
  //  }
}