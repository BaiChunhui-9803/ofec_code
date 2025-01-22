/********* Begin Register Information **********
{
    "description": "niche center distinguish-based differential evolution",
    "identifier": "NCD_DE",
    "name": "NCD-DE",
    "tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/**********************************************
@article{jiang2023optimizing,
  volume = {53},
  number = {4},
  journal = {IEEE Transactions on Cybernetics},
  pages = {2544--2557},
  year = {2023},
  author = {Yi Jiang and Zhihui Zhan and Kay Chen Tan and Jun Zhang},
  title = {Optimizing niche center for multimodal optimization problems}
}
***********************************************/

#ifndef OFEC_NCD_DE_H
#define OFEC_NCD_DE_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/differential_evolution/population.h"

namespace ofec {
    class NicheCenterDistinguish : public ProblemVariableVector<int> {
        OFEC_CONCRETE_INSTANCE(NicheCenterDistinguish)
    protected:
        const PopulationDE<> *m_pop = nullptr;
        Problem *m_original_problem = nullptr;
        std::vector<std::vector<Real>> m_dis;
  
        void addInputParameters() {}
        void initialize_(Environment *env) override;
        void evaluate(const VariableBase &vars, std::vector<Real> &objs, std::vector<Real> &cons) const override;

    public:
        void setData(const PopulationDE<> &pop, Problem *original_problem);
        void initializeVariables(VariableBase &x, Random *rnd) const override;
    };

    class InternalGA : virtual public Algorithm {
        OFEC_CONCRETE_INSTANCE(InternalGA)
    protected:
        size_t m_pop_size;
        Real m_crossover_rate, m_mutation_rate;
        size_t m_iteration_epochs;
        Population<Solution<VariableVector<int>>> m_pop;
        std::vector<size_t> m_rand_seq; // Random sequence of the population for tournamentSelection()

        void addInputParameters() {}
        void initialize_(Environment *env) override;
        void run_(Environment *env) override;
        size_t tournamentSelection(const Population<Solution<VariableVector<int>>> &pop, Environment *env);
        void crossover(const Solution<VariableVector<int>> &p1, const Solution<VariableVector<int>> &p2,
            Solution<VariableVector<int>> &o1, Solution<VariableVector<int>> &o2, Environment *env);
        void mutate(Solution<VariableVector<int>> &s);

    public:
        const VariableVector<int>& bestChromosome(Environment *env);
    };

	class NCD_DE : virtual public Algorithm {
        OFEC_CONCRETE_INSTANCE(NCD_DE)
	protected:
        size_t m_pop_size, m_lambda;
        Real m_scaling_factor, m_crossover_rate;
        PopulationDE<> m_pop;

        void addInputParameters();
        void initialize_(Environment *env) override;
		void run_(Environment *env) override;
        size_t nearestInd(const Solution<> &s, Environment *env);
	};
}

#endif // ! OFEC_NCD_DE_H
