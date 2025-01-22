#ifndef OFEC_GA_ONEMAX_H
#define OFEC_GA_ONEMAX_H

#include "geneticAlgorithm.hpp"
#include "../../../../../core/algorithm/algorithm.h"
#include "../../../../problem/combination/one_max/one_max.h"
namespace ofec{
 

        class GeneticAlgorithmOneMax : protected modularGA::GeneticAlgorithm, public Algorithm {
        protected:

            std::vector< OneMax::solutionType> m_totalSols;
            std::unique_ptr<OneMax::solutionType> m_sol;



            virtual double Evaluate(std::vector<int>& x) override {
                using namespace std;
                double result;
                m_sol->variable().vect() = x;
                m_sol->evaluate(m_problem.get(), this);

                ofec::Real pos = (m_problem->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
                m_sol->setFitness(pos * m_sol->objective()[0]);
              //  m_sol->transferObjToFit(m_problem.get());


                result = m_sol->fitness();

                // if (this->csv_logger_ != nullptr) {
                //   this->csv_logger_->do_log(this->problem_->loggerInfo()); /// TODO: we assume only using PBO suite now.
                // }
                ++this->evaluation_;

                if (result > this->best_found_fitness_) {
                    this->best_found_fitness_ = result;
                    this->best_individual_ = x;
                }



                // This is only used to calculate ERT
                if (this->ERT_flag_) {
                    /*if (Opt == optimizationType::MAXIMIZATION) */ {
                        if (!this->hitting_flag_ && this->best_found_fitness_ >= this->hiting_target_) {
                            this->hitting_time_.push_back(this->evaluation_);
                            this->hitting_flag_ = true;
                        }
                    }
                    //else {
                    //    if (!this->hitting_flag_ && this->best_found_fitness_ <= this->hiting_target_) {
                    //        this->hitting_time_.push_back(this->evaluation_);
                    //        this->hitting_flag_ = true;
                    //    }
                    //}
                }
                return result;
            }

            virtual void evolve();
            virtual void Initialization() {
                int n = m_problem->numVariables();

                for (int i = 0; i != this->mu_; ++i) {
                    std::vector<int> tmp(n, 0);
                    for (int j = 0; j != n; ++j) {
                        if (m_random->uniform.next() < 0.5) {
                            tmp[j] = 1;
                        }
                    }
                    this->parents_population_.push_back(tmp);
                    this->parents_fitness_.push_back(this->Evaluate(tmp));
                }
            }
        public:

            void run_() override;
            void initialize_() override;
#ifdef OFEC_DEMO
            void updateBuffer();
#endif
        };
    
}


#endif //  GA_ONEMAX_H
