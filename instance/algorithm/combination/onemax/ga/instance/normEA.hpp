/********* Begin Register Information **********
{
    "name": "NormEA-OneMax",
    "identifier": "NormEA_OneMax",
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
#include "../../../../../../utility/typevar/typevar.h"
namespace ofec
{


      class NormEA_OneMax : public GeneticAlgorithmOneMax
      {
          double rand;
          std::vector<int> best_ind, x;
          double r_n ;
          double best_f_x, f_x;
          int flip_n;
          double init_r_n_ = 1.0;
      public:
          static void setParameter(ParamMap& param, int lambda, double init_r = 1.0)
        {


              param["lambda"] = lambda;
              param["init_r"] = init_r;

        //    init_r_n_ = init_r;
          //  
          //this->set_mu(1);
          //this->set_lambda(lambda);
          //this->set_mutation_operator("NORMALSAMPLE");
          //this->set_r_n(init_r_n_);
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
            auto& v = *m_param;
            lambda = v.get<int>("lambda");
            init_r_n_ = v.get<double>("init_r");


            this->set_mu(1);
            this->set_lambda(lambda);
            this->set_mutation_operator("NORMALSAMPLE");
            this->set_r_n(init_r_n_);


            this->Preparation(m_problem->numVariables());
            this->Initialization();
            r_n = this->init_r_n_;
        }

        //void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
        //{
        //  string algorithm_name = "(1+" + to_string(this->get_lambda()) + ")-normEA>0";
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


        virtual void evolve() override{
  

            this->update_generation();

            best_f_x = std::numeric_limits<int>::min();
            this->set_r_n(r_n);
            this->set_sigma_n(sqrt(static_cast<double>(r_n) * (1.0 - r_n / static_cast<double>(m_problem->numVariables()))));
            for (size_t i = 0; i < this->get_lambda(); ++i)
            {
                x = this->get_parents_population()[0];
                flip_n = this->DoMutation(x, m_random.get());

                f_x = this->Evaluate(x);
                if (f_x > best_f_x)
                {
                    best_f_x = f_x;
                    r_n = flip_n;
                    best_ind = x;
                }

                //if (this->Termination())
                //    break;
            }

            //if (this->Termination())
            //    break;

            if (best_f_x >= this->get_parents_fitness()[0])
            {
                this->set_parents_population(best_ind, 0);
                this->set_parents_fitness(best_f_x, 0);
            }
        }
        //void DoGeneticAlgorithm()
        //{
        //    using namespace std;
        //  double rand;
        //  vector<int> best_ind, x;
        //  double r_n = this->init_r_n_;
        //  double best_f_x, f_x;
        //  int flip_n;
        //  this->Preparation();
        //  this->Initialization();

        //  while (!this->Termination())
        //  {
        //    this->update_generation();

        //    best_f_x = std::numeric_limits<int>::min();
        //    this->set_r_n(r_n);
        //    this->set_sigma_n(sqrt(static_cast<double>(r_n) * (1.0 - r_n/static_cast<double>(this->get_dimension()))));
        //    for (size_t i = 0; i < this->get_lambda(); ++i)
        //    {
        //      x = this->get_parents_population()[0];
        //      flip_n = this->DoMutation(x);
        //      
        //      f_x = this->Evaluate(x);
        //      if (f_x > best_f_x)
        //      {
        //        best_f_x = f_x;
        //        r_n = flip_n;
        //        best_ind = x;
        //      }

        //      if (this->Termination())
        //        break;
        //    }

        //    if (this->Termination())
        //      break;

        //    if (best_f_x >= this->get_parents_fitness()[0])
        //    {
        //      this->set_parents_population(best_ind, 0);
        //      this->set_parents_fitness(best_f_x, 0);
        //    }

        //  }
        //}

        protected:

        //std::vector<int> best_ind;
        ////double r_n = this->init_r_n_;
        //double best_f_x;

      };
    
}