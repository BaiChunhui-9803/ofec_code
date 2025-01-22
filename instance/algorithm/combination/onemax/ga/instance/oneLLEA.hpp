/********* Begin Register Information **********
{
    "name": "oneLambdaLambdaEA-OneMax",
    "identifier": "oneLambdaLambdaEA_OneMax",
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

 // namespace ga
 // {

      class oneLambdaLambdaEA_OneMax : public GeneticAlgorithmOneMax
      {
      protected:
          double a = std::pow(1.5, 0.25); /// < parameter for adjusting lambda
          double b = 2.0 / 3.0;      /// < parameter for adjusting lambda
          double best_f = std::numeric_limits<double>::lowest();


          double rand, best_mutation_f;
          std::vector<int> mutation_offspring, best_ind, x;
          bool update_lambda_flag;
          int lambda, dimension;
          double mutation_rate, p_u;
          // int mutation_strength;



      public:
          static void setParameter(ParamMap& param, int lambda)
        {
              param["lambda"] = lambda;
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

              //{
              //    VisualParamterRealSlider curPara;
              //    curPara.m_key = "init_r";
              //    curPara.setRange(1, 30, 1, 0.1);
              //    param.emplace_back(new VisualParamterRealSlider(curPara));
              //}
          }


          void initialize_() override {
              GeneticAlgorithmOneMax::initialize_();
              int lambda;
              auto& v = *m_param;
              lambda = v.get<int>("lambda");


              this->set_mu(1);
              this->set_lambda(lambda);
              this->set_mutation_operator("BINOMIALSAMPLE");
              this->set_crossover_operator("UNIFORMCROSSOVER");
              this->set_selection_operator("BESTPLUS");


              a = std::pow(1.5, 0.25); /// < parameter for adjusting lambda
              b = 2.0 / 3.0;      /// < parameter for adjusting lambda
              best_f = std::numeric_limits<double>::lowest();
              this->Preparation(m_problem->numVariables());

              this->Initialization();


              lambda = this->get_lambda();
              dimension = m_problem->numVariables();
          }
        //void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective> > suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
        //{
        //  string algorithm_name = "(1+(" + to_string(this->get_lambda()) + "," + to_string(this->get_lambda()) + "))>_0 EA";
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


 

         // while (!this->Termination())
          {

            this->update_generation();
            update_lambda_flag = false;
            best_mutation_f = std::numeric_limits<double>::lowest();
            this->clear_offspring_population();
            this->clear_offspring_fitness();

            mutation_rate = static_cast<double>(lambda) / static_cast<double>(dimension);
            this->set_p_u(1.0 / static_cast<double>(lambda));

            /**Mutation stage
             */
           auto mutation_strength = this->SampleConditionalBinomial(mutation_rate, dimension, m_random.get());
            for (size_t i = 0; i < this->get_lambda(); ++i)
            {
              x = this->get_parents_population()[0];

              this->Flip(x, mutation_strength,m_random.get());
              this->add_offspring_fitness(this->Evaluate(x));

              if (this->get_offspring_fitness()[i] >= best_f)
              {
                best_ind = x;
                best_f = this->get_offspring_fitness()[i];
                if (this->get_offspring_fitness()[i] > best_f)
                {
                  update_lambda_flag = true;
                }
              }

              if (this->get_offspring_fitness()[i] > best_mutation_f)
              {
                best_mutation_f = this->get_offspring_fitness()[i];
                mutation_offspring = x;
              }

          //    if (this->Termination())
            //    break;
            }

            /** Crossover Stage
             */
            for (size_t i = 0; i < this->get_lambda(); ++i)
            {

              this->DoCrossover(x, this->get_parents_population()[0], mutation_offspring,m_random.get());
              if (this->c_flipped_index.size() == 0)
              {
                this->set_offspring_fitness(this->get_parents_fitness()[0], i);
              }
              else if (x == mutation_offspring)
              {
                this->set_offspring_fitness(best_mutation_f, i);
              }
              else
              {
                this->set_offspring_fitness(this->Evaluate(x), i);
              }

              if (this->get_offspring_fitness()[i] >= best_f)
              {
                best_ind = x;
                best_f = this->get_offspring_fitness()[i];
                if (this->get_offspring_fitness()[i] > best_f)
                {
                  update_lambda_flag = true;
                }
              }

              //if (this->Termination())
              //  break;
            }

         //   if (this->Termination())
         //     break;

            this->set_parents_population(best_ind, 0);
            this->set_parents_fitness(best_f, 0);
            if (update_lambda_flag == true)
            {
              lambda = (lambda * b) > 1 ? (lambda * b) : 1;
            }
            else
            {
              lambda = (lambda * a) < dimension ? (lambda * a) : dimension;
            }
          }
        }
      };
    
//  }
}