/********* Begin Register Information **********
{
	"name": "DCNSGAII-DE",
	"identifier": "DCNSGAII_DE",
	"problem tags": [ "ConOP", "COP", "MOP", "SOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_DCNSGAII_DE_H
#define OFEC_DCNSGAII_DE_H

#include "../../../../.././../core/algorithm/population.h"
#include "../../../../template/classic/de/Solution.h"
#include "../../../../template/framework/dcmoea/dcmoea.h"

namespace ofec {
	class DCNSGAII_DE_pop : public Population<DCMOEA_ind<IndDE>>, DCMOEA<DCMOEA_ind<IndDE>> {
	public:
		using typename Population<DCMOEA_ind<IndDE>>::IndividualType;
	public:
		DCNSGAII_DE_pop(size_t size_pop, Problem *pro);
		int evolve(Problem *pro, Algorithm *alg, Random *rnd);
		void initialize(Problem *pro, Algorithm *alg, Random *rnd);
		void sort(Problem *pro);
	protected:
		void select_next_parent_population(std::vector<IndividualType*>& pop, Problem *pro);
		void reduce_radius(Problem *pro);
		void caculate_nichecount(std::vector<IndividualType*>& pop, Problem *pro);
		Dominance e_Pareto_compare(IndividualType* const& s1, IndividualType* const& s2, Problem *pro);
		Dominance Pareto_compare(IndividualType* const& s1, IndividualType* const& s2, Problem *pro);
	protected:
		Real m_r = 0, m_R, m_crossover_rate, m_scaling_factor;
		std::vector<std::unique_ptr<IndividualType>> m_offspring;
		std::vector<std::unique_ptr<IndividualType>> m_temp_pop;
		RecombineDE m_recombine_strategy;
		std::vector<int> m_rand_seq;
		bool m_flag = false;
	};

	class DCNSGAII_DE : public Algorithm {
	protected:
		std::unique_ptr<DCNSGAII_DE_pop> m_pop;
		size_t m_pop_size;
	public:
		void initialize_();
		void run_();
		void initPop();
		void record();
	};
}



#endif // !OFEC_DCNSGAII_DE_H
