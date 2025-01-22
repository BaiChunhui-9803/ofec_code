#include "dynde_subpopulation.h"

namespace ofec {
	SubPopDynDE::SubPopDynDE(size_t size, Problem *pro) : 
		PopDE<IndDynDE>(size, pro), 
		m_num_normal(5), 
		m_num_brownian(5), 
		m_num_quantum(5),
		m_r_cloud(1.0), 
		m_sigma(0.2) 
	{
		// set m_Rcloud to shift lenght, default for shift length =1.0 if not known 
		//default configuration (N,N+) or (N,Nq)= (5,5)
		assign_type();
		m_mutation_strategy =  MutationDE::best_2;
		//update_archive(*m_individuals[0]);
	}

	void SubPopDynDE::assign_type() {
		// for the version with the suggested settings 
		for (size_t i = 0; i < size(); ++i) {
			if (i < m_num_normal) m_individuals[i]->m_type = IndDynDE::Solution_type::TYPE_DE;
			else  m_individuals[i]->m_type = IndDynDE::Solution_type::TYPE_BROWNIAN;
		}
	}


	int SubPopDynDE::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		if (size() < 1)
			return kNormalEval;
		int tag = kNormalEval;
		updateBest(pro);
		for (size_t i = 0; i < size(); ++i) {
			if (m_individuals[i]->m_type == IndDynDE::Solution_type::TYPE_DE || m_individuals[i]->m_type == IndDynDE::Solution_type::TYPE_ENTROPY_DE) {
				mutate(i, rnd, pro);
				m_individuals[i]->recombine(m_crossover_rate, m_recombine_strategy, rnd, pro);
			}
			else if (m_individuals[i]->m_type == IndDynDE::Solution_type::TYPE_BROWNIAN) {
				tag = m_individuals[i]->brownian(*m_best.front(), m_sigma, pro, alg, rnd);

			}
			else if (m_individuals[i]->m_type == IndDynDE::Solution_type::TYPE_QUANTUM) {
				tag = m_individuals[i]->quantum(*m_best.front(), m_r_cloud, pro, alg, rnd);
			}
			if (!(tag & kNormalEval))  return tag;
		}
		for (int i = 0; i < size(); i++) {
			if (m_individuals[i]->m_type == IndDynDE::Solution_type::TYPE_DE || m_individuals[i]->m_type == IndDynDE::Solution_type::TYPE_ENTROPY_DE) {
				tag = m_individuals[i]->select(pro, alg);
				if (!(tag & kNormalEval)) return tag;
			}
			if (m_individuals[i]->m_type == IndDynDE::Solution_type::TYPE_ENTROPY_DE) {
				// add entropy if necessary 
				tag = m_individuals[i]->entropy(m_sigma, pro, alg, rnd);
				if (!(tag & kNormalEval)) return tag;
			}
		}
		if (tag & kNormalEval) {
			m_iteration++;
		}
		return tag;
	}
}
