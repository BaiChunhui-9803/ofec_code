#ifndef OFEC_GA_POP_H
#define OFEC_GA_POP_H

#include "../../../../../core/algorithm/population.h"
#include <numeric>
#include "../../../../../utility/functional.h"

namespace ofec {
	template<typename TInd>
	class PopGA : public Population<TInd> {
	private:
		std::vector<size_t> m_rand_seq; // Random sequence of the population for tournamentSelection()
	protected:
		Real m_cr = 0.9; // crossover probability
		Real m_mr = 0.1; // mutation probability

	public:
		PopGA() = default;
		template<typename ... Args>
		PopGA(size_t size_pop, Environment *env, Args&& ...args) :
			Population<TInd>(size_pop, env, std::forward<Args>(args)...) {}
		void setRate(Real cr, Real mr) { m_cr = cr; m_mr = mr; }
		void setCrossRate(Real cr) { m_cr = cr; }
		void setMutateRate(Real mr) { m_mr = mr; }
		virtual void crossover(size_t idx_parent1, size_t idx_parent2, TInd &child1, TInd &child2, Environment *env, Random* rnd) = 0;
		virtual void mutate(TInd &ind, Environment *env, Random* rnd) = 0;
		size_t tournamentSelection(Environment *env, Random* rnd, size_t tournament_size = 2);
		void reproduction(Population<TInd> &pop_and_off, Environment *env, Random* rnd);
	};

	template<typename TInd>
	size_t PopGA<TInd>::tournamentSelection(Environment *env, Random *rnd, size_t tournament_size) {
		if (m_rand_seq.size() != this->m_individuals.size()) {
			m_rand_seq.resize(this->m_individuals.size());
			std::iota(m_rand_seq.begin(), m_rand_seq.end(), 0);
		}
		rnd->uniform.shuffle(m_rand_seq.begin(), m_rand_seq.end());
		size_t idx_best = m_rand_seq[0];
		for (size_t k = 1; k < tournament_size; ++k)
			if (dominate(*this->m_individuals[m_rand_seq[k]], *this->m_individuals[idx_best], env->problem()->optimizeMode())) {
				idx_best = m_rand_seq[k];
			}
		return idx_best;
	}

	template<typename TInd>
	void PopGA<TInd>::reproduction(Population<TInd> &pop_and_off, Environment *env, Random *rnd) {
		std::vector<size_t> p(2);
		for (size_t i = 0; i < this->m_individuals.size(); i += 2) {
			p[0] = tournamentSelection(env, rnd);
			do { p[1] = tournamentSelection(env, rnd); } while (p[1] == p[0]);
			crossover(p[0], p[1], pop_and_off[i], pop_and_off[i + 1], env, rnd);
			mutate(pop_and_off[i], env, rnd);
			mutate(pop_and_off[i + 1], env, rnd);
		}

		for (size_t i = 0; i < this->m_individuals.size(); ++i) {
			pop_and_off[i + this->m_individuals.size()] = *this->m_individuals[i];
		}

	}


}

#endif // !OFEC_GA_POP_H

