#include "ga_onemax.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

void ofec::GeneticAlgorithmOneMax::evolve()
{
    this->offspring_population_.clear();
    this->offspring_fitness_.clear();
    for (size_t i = 0; i < this->lambda_; ++i) {
        this->SelectTwoParents(m_random.get());
        this->offspring_population_.push_back(this->parents_population_[this->selected_parents_[0]]);

        auto rand = m_random->uniform.next();
        this->c_flipped_index.clear();
        this->m_flipped_index.clear();
        if (rand < this->crossover_probability_) {
            this->DoCrossover(this->offspring_population_[i], this->parents_population_[this->selected_parents_[0]], this->parents_population_[this->selected_parents_[1]],m_random.get());
        }

        if (this->crossover_mutation_r_) {
            this->DoMutation(this->offspring_population_[i],m_random.get());
        }
        else if (rand >= this->crossover_probability_) {
            this->DoMutation(this->offspring_population_[i],m_random.get());
        }

        if (this->c_flipped_index == this->m_flipped_index) { /// If the flipping indexes of crossover and mutation are identical, the individual remains the same.
            this->offspring_fitness_.push_back(this->parents_fitness_[this->selected_parents_[0]]);
        }
        else if (rand < this->crossover_probability_ && this->offspring_population_[i] == this->parents_population_[this->selected_parents_[1]]) { /// If the offspring is identical with the second parent.
       /// TODO: Do something to save time for this comparison.
            this->offspring_fitness_.push_back(this->parents_fitness_[this->selected_parents_[1]]);
        }
        else { /// otherwise evaluate.
            this->offspring_fitness_.push_back(this->Evaluate(this->offspring_population_[i]));
        }

        //    if (this->Termination()) break;
    }

    //if (this->Termination()) break;

    this->DoSelection(this->parents_population_, this->parents_fitness_, this->offspring_population_, this->offspring_fitness_,m_random.get());
    this->AdaptiveStrategy();
}

void ofec::GeneticAlgorithmOneMax::run_()
{
    while (!terminating()) {
        evolve();

        if (m_best_worst_monitored) {
            std::unique_ptr<ofec::SolBase> bestSol(m_problem->createSolution());
            bestSol->setFitness(-std::numeric_limits<Real>::max());
            for (auto& it : parents_population_) {
                m_sol->variable() = it;
                m_sol->evaluate(m_problem.get(),nullptr);

                ofec::Real pos = (m_problem->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
                m_sol->setFitness(pos * m_sol->objective()[0]);
             //   m_sol->transferObjToFit(m_problem.get());
                if (bestSol->fitness() < m_sol->fitness()) {
                    bestSol.reset(m_problem->createSolution(*m_sol));
                }

               // m_totalSols.push_back(*m_sol);

            }
            m_best_so_far.emplace_back(std::move(bestSol));
        }

#ifdef OFEC_DEMO
        updateBuffer();
#endif
    }
}

void ofec::GeneticAlgorithmOneMax::initialize_()
{
    using namespace modularGA;
    Algorithm::initialize_();
    

    mu_ = DEFAULT_MU_;
    lambda_ = (DEFAULT_LAMBDA_);
    crossover_probability_ = (DEFAULT_CROSSOVER_PROBABLITY_);
    crossover_mutation_r_ = (DEFAULT_CROSSOVER_MUTATION_RELATION_);
    evaluation_ = 0;
    generation_ = 0;
    evluation_budget_ = DEFAULT_EVALUATION_BUDGET_;
    generation_budget_ = DEFAULT_GENERATION_BUDGET_;
    optimum_ = std::numeric_limits<double>::max();/// < TODO: now we assume doing maximization.

    m_sol.reset(new OneMax::solutionType(m_problem->numObjectives(), m_problem->numConstraints(), m_problem->numVariables()));

    Preparation(m_problem->numVariables());

    Initialization();

    if (m_best_worst_monitored) {
        m_best_so_far.clear();
        
    }

   

}
#ifdef OFEC_DEMO
void ofec::GeneticAlgorithmOneMax::updateBuffer()
{
    
    m_solution.clear();
    m_solution.resize(1);
    m_totalSols.clear();
    for (auto& it : parents_population_) {
        m_sol->variable() = it;
        m_totalSols.push_back(*m_sol);
        
    }

    for (auto& it : offspring_population_) {
        m_sol->variable() = it;
        m_totalSols.push_back(*m_sol);

    }

    for (size_t i = 0; i < m_totalSols.size(); ++i)
        m_solution[0].push_back(&m_totalSols.at(i));
    ofec_demo::g_buffer->appendAlgBuffer(this);
}
#endif
