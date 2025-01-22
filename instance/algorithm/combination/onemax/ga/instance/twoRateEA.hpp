/********* Begin Register Information **********
{
    "name": "TwoRateEA-OneMax",
    "identifier": "TwoRateEA_OneMax",
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
namespace ofec
{

//  namespace ga
 // {
      class TwoRateEA_OneMax : public GeneticAlgorithmOneMax
      {
          double r = 0;
          std::vector<int> best_ind, x;

          double f_x, best_f_x;
          size_t best_x_index;
          //= this->init_r_;
      public:
          static void setParameter(ParamMap& param,int lambda, double init_r = 2.0)
        {
              param["lambda"] = lambda;
              param["init_r"] = init_r;
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
                  curPara.m_key = "init_r";
                  curPara.setRange(1, 30, 1, 0.1);
                  param.emplace_back(new VisualParamterRealSlider(curPara));
              }
          }

        void initialize_() override {
            GeneticAlgorithmOneMax::initialize_();


            int lambda;
            double init_r;
            auto& v = *m_param;
            lambda = v.get<int>("lambda");
            init_r = v.get<double>("init_r");
            this->set_mu(1);
            this->set_lambda(lambda);
            this->set_crossover_mutation_r("OR");
            this->set_crossover_probability(0);
            this->set_mutation_operator("BINOMIALSAMPLE");
            this->set_selection_operator("BESTPLUS");
            this->init_r_ = init_r;


            this->Preparation(m_problem->numVariables());
            this->Initialization();

            r = this->init_r_;
        }

        //void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
        //{
        //  string algorithm_name = "(1+" + to_string(this->get_lambda()) + ")-2rate-EA>0";
        //  std::shared_ptr<ioh::logger::Analyzer > logger(new ioh::logger::Analyzer(
        //                                              {ioh::trigger::on_improvement}, // trigger when the objective value improves
        //                                              {},                   // no additional properties 
        //                                              folder_path,        // path to store data
        //                                              folder_name,               // name of the folder in path, which will be newly created
        //                                              algorithm_name,                     // name of the algoritm 
        //                                              algorithm_name,                   // additional info about the algorithm              
        //                                              false            // where to store x positions in the data files 
        //                                            ));                                                  
        //  this->AssignLogger(logger);

        //  this->set_evaluation_budget(eval_budget);
        //  this->set_generation_budget(gene_budget);
        //  this->set_independent_runs(independent_runs);
        //  this->SetSeed(rand_seed);

        //  this->run_N(suite);
        //}

        virtual void evolve() override
        {
          //double rand;



       //   while (!this->Termination())
          {
            this->update_generation();

            best_f_x = std::numeric_limits<int>::min();
            best_x_index = 0;
            for (size_t i = 0; i < this->get_lambda(); ++i)
            {

              if (i < static_cast<size_t>(floor(this->get_lambda() / 2.0)))
              {
                x = this->get_parents_population()[0];
                this->set_mutation_rate(r / 2.0 / static_cast<double>(m_problem->numVariables()));
                this->DoMutation(x,m_random.get());
              }
              else
              {
                x = this->get_parents_population()[0];
                this->set_mutation_rate(r * 2.0 / static_cast<double>(m_problem->numVariables()));
                this->DoMutation(x, m_random.get());
              }

              f_x = this->Evaluate(x);
              if (f_x > best_f_x)
              {
                best_f_x = f_x;
                best_x_index = i;
                best_ind = x;
              }

  /*            if (this->Termination())
                break;*/
            }

     /*       if (this->Termination())
              break;*/

            if (best_f_x >= this->get_parents_fitness()[0])
            {
              this->set_parents_population(best_ind, 0);
              this->set_parents_fitness(best_f_x, 0);
            }

            if (m_random->uniform.next() < 0.5)
            {
              if (best_x_index < static_cast<size_t>(floor(this->get_lambda() / 2.0)))
              {
                r = r / 2.0;
              }
              else
              {
                r = r * 2.0;
              }
            }
            else
            {
              if (m_random->uniform.next() < 0.5)
              {
                r = r / 2.0;
              }
              else
              {
                r = r * 2.0;
              }
            }
            r = r < 2.0 ? 2.0 : r;
            r = r > (double(m_problem->numVariables()) / 4.0) ? (double(m_problem->numVariables()) / 4.0) : r;
          }
        }

        double init_r_ = 2.0;
      };
      //   }
}