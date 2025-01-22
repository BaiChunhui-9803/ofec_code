#ifndef OFEC_JADE_POP_H
#define OFEC_JADE_POP_H

#include "../../../../template/classic/differential_evolution/population.h"

namespace ofec {
	template<typename TInd = IndividualDE>
	class PopJADE : public PopulationDE<TInd> {
	public:
		PopJADE();
		PopJADE(size_t size_pop, Environment *env);
		void resize(size_t size_pop, Environment *env) override;
		void initialize(Environment *env, Random *rnd) override;
		int evolve(Environment *env, Random *rnd) override;
		Real F(size_t i) const { return mv_F[i]; }
		Real CR(size_t i) const { return mv_CR[i]; }

	protected:
		void selectTrial(size_t base, Random *rnd);
		void updateF(Random *rnd);
		void updateCR(Random *rnd);

	protected:
		///////////////////////algorithm parameters/////////////////////////////
		Real m_p;
		Real m_c;
		int m_archive_flag;
		std::vector<Solution<>> m_archive;
		std::vector<Solution<>*> m_candi;
		std::vector<Real> m_pcent_best;
		std::vector<int> m_pcent_best_index;
		std::vector<Real> mv_F, mv_CR;
		Real m_mu_F, m_mu_CR;
		///////////////////////////////////////////////////////////////////////////
	};

	template<typename TInd>
	PopJADE<TInd>::PopJADE() :
		PopulationDE<TInd>(),
		m_candi(3, nullptr),
		m_archive_flag(1),
		m_p(0.2),
		m_c(0.1),
		m_mu_F(0.5),
		m_mu_CR(0.5) {}

	template<typename TInd>
	PopJADE<TInd>::PopJADE(size_t size_pop, Environment *env) :
		PopulationDE<TInd>(size_pop, env),
		m_candi(3, nullptr),
		m_pcent_best(size_pop),
		m_pcent_best_index(size_pop),
		mv_F(size_pop, 0),
		mv_CR(size_pop, 0),
		m_archive_flag(1),
		m_p(0.2),
		m_c(0.1),
		m_mu_F(0.5),
		m_mu_CR(0.5) {}

	template<typename TInd>
	void PopJADE<TInd>::resize(size_t size_pop, Environment *env) {
		PopulationDE<TInd>::resize(size_pop, env);
		m_pcent_best.resize(size_pop);
		m_pcent_best_index.resize(size_pop);
		mv_F.resize(size_pop);
		mv_CR.resize(size_pop);
	}

	template<typename TInd>
	void PopJADE<TInd>::initialize(Environment *env, Random *rnd) {
		Population<TInd>::initialize(env, rnd);
		m_mu_CR = 0.5;
		m_mu_F = 0.5;
		m_archive.clear();
		m_archive_flag = 1;
	}

	template<typename TInd>
	void PopJADE<TInd>::selectTrial(size_t base, Random *rnd) {
		std::vector<int> candidate;
		for (size_t i = 0; i < this->size(); i++) {
			if (base != i) candidate.push_back(i);
		}
		size_t idx;
		do {
			idx = m_pcent_best_index[rnd->uniform.nextNonStd<size_t>(0, this->size() * m_p)];
		} while (idx == base);
		m_candi[0] = this->m_individuals[idx].get();
		const auto it = std::find(candidate.begin(), candidate.end(), idx);
		candidate.erase(it);
		idx = rnd->uniform.nextNonStd<int>(0, (int)candidate.size());
		m_candi[1] = this->m_individuals[candidate[idx]].get();
		candidate.erase(candidate.begin() + idx);
		idx = rnd->uniform.nextNonStd<int>(0, (int)(candidate.size() + m_archive.size()));
		if (idx >= candidate.size()) {
			m_candi[2] = &m_archive[idx - candidate.size()];
		}
		else {
			m_candi[2] = this->m_individuals[candidate[idx]].get();
		}
	}

	template<typename TInd>
	int PopJADE<TInd>::evolve(Environment *env, Random *rnd) {
		updateCR(rnd);
		updateF(rnd);
		int tag = kNormalEval;
		for (size_t i = 0; i < this->size(); ++i) {
			m_pcent_best[i] = this->m_individuals[i]->objective()[0];
		}
		mergeSort(m_pcent_best, (int)this->size(), m_pcent_best_index, env->problem()->optimizeMode(0) == OptimizeMode::kMinimize);
		std::vector<int> candidate(3);
		for (size_t i = 0; i < this->size(); ++i) {
			selectTrial(i, rnd);
			this->m_individuals[i]->mutate(mv_F[i], this->m_individuals[i].get(), m_candi[0], 
				this->m_individuals[i].get(), env, m_candi[1], m_candi[2]);
			this->m_individuals[i]->recombine(mv_CR[i], this->m_recombine_strategy, rnd, env);
		}
		for (size_t i = 0; i < this->size(); i++) {
			tag = this->m_individuals[i]->select(env);
			if (!this->m_individuals[i]->isImproved()) {
				m_archive.push_back(*this->m_individuals[i]);
			}
			while (m_archive.size() > this->size()) {
				int ridx = rnd->uniform.nextNonStd<int>(0, m_archive.size());
				m_archive.erase(m_archive.begin() + ridx);
			}
			if (!(tag & kNormalEval)) {
				return tag;
			}
		}
		if (tag & kNormalEval) {
			this->m_iteration++;
		}
		return tag;
	}

	template<typename TInd>
	void PopJADE<TInd>::updateF(Random *rnd) {
		if (this->m_iteration > 0) {
			Real mean = 0, mean2 = 0;
			int cnt = 0;
			for (size_t i = 0; i < this->size(); i++) {
				if (this->m_individuals[i]->isImproved()) {
					cnt++;
					mean += mv_F[i];
					mean2 += mv_F[i] * mv_F[i];
				}
			}
			if (cnt > 0) {
				mean = mean2 / mean;
				m_mu_F = (1 - m_c) * m_mu_F + m_c * mean;
			}
			else { // Junchen Wang's personal opinion as the case is not mentioned in the paper 
				m_mu_F = rnd->uniform.next();
			}
		}
		for (size_t i = 0; i < this->size(); i++) {
			do {
				mv_F[i] = rnd->cauchy.nextNonStd(m_mu_F, 0.1);
			} while (mv_F[i] <= 0);
			if (mv_F[i] > 1) mv_F[i] = 1;
		}
	}

	template<typename TInd>
	void PopJADE<TInd>::updateCR(Random *rnd) {
		if (this->m_iteration > 0) {
			Real mean = 0;
			int cnt = 0;
			for (size_t i = 0; i < this->size(); i++) {
				if (this->m_individuals[i]->isImproved()) {
					cnt++;
					mean += mv_CR[i];
				}
			}
			if (cnt > 0) {
				mean /= cnt;
				m_mu_CR = (1 - m_c) * m_mu_CR + m_c * mean;
			}
			else { // Junchen Wang's personal opinion as the case is not mentioned in the paper 
				m_mu_CR = rnd->uniform.next();
			}
		}
		for (size_t i = 0; i < this->size(); i++) {
			mv_CR[i] = rnd->normal.nextNonStd(m_mu_CR, 0.1);
			if (mv_CR[i] < 0) mv_CR[i] = 0;
			else if (mv_CR[i] > 1) mv_CR[i] = 1;
		}
	}
}
#endif // OFEC_JADE_POP_H
