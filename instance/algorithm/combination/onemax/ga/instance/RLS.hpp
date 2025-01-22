/********* Begin Register Information **********
{
    "name": "RLS-OneMax",
    "identifier": "RLS_OneMax",
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
//  {
      class RLS_OneMax : public GeneticAlgorithmOneMax
      {

      public:
          static void setParameter(ParamMap& param)
        {
          //this->set_mu(1);
          //this->set_lambda(1);
          //this->set_crossover_mutation_r("IND");
          //this->set_crossover_probability(0);
          //this->set_mutation_operator("STATICSAMPLE");
          //this->set_l(1);
          //this->set_selection_operator("BESTPLUS");
        }

          void initialize_() override {
              GeneticAlgorithmOneMax::initialize_();
              this->set_mu(1);
              this->set_lambda(1);
              this->set_crossover_mutation_r("IND");
              this->set_crossover_probability(0);
              this->set_mutation_operator("STATICSAMPLE");
              this->set_l(1);
              this->set_selection_operator("BESTPLUS");


          }

        //void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective> > suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
        //{
        //  string algorithm_name = "RLS";
        //   std::shared_ptr<ioh::logger::Analyzer > logger(new ioh::logger::Analyzer(
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
      };
  //  }
}