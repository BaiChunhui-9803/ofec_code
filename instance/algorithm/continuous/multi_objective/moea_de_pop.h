/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*  See the details of MOEA/D-DE in the following paper
*  H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization
*  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
*  University of Essex, 2009
*************************************************************************/
// Created: 30 December 2014
// Last modified: 18 July 2019 by Xiaofang Wu (email:wuxiaofang@cug.edu.cn)

#ifndef MOEA_DE_POP_H
#define MOEA_DE_POP_H
#include "../../template/classic/de/individual.h"
#include "../../template/classic/de/population.h"
#include <numeric>

namespace ofec {
	template<typename TInd = IndDE>
	class PopMODE : public PopDE<TInd> {
	protected:
		Real m_mr;
		Real m_meta;
		std::vector<size_t> m_rand_seq; // Random sequence of the population

	public:
		PopMODE(size_t size_pop, Problem *pro);
		//PopMODE<TInd>& operator=(const PopMODE& rhs);
		void setParamMODE(Real mr, Real meta) { m_mr = mr; m_meta = meta; }
		void crossMutate(const std::vector<size_t>& index, TInd& child, Problem *pro, Random *rnd);
		size_t tournamentSelection(Problem *pro, Random *rnd, size_t tournament_size = 2);

	protected:
		void crossover(const TInd& parent1, const TInd& parent2, const TInd& parent3, TInd& child, Problem *pro, Random *rnd);
		void mutate(TInd& child, Problem *pro, Random *rnd);
	};

	template<typename TInd>
	PopMODE<TInd>::PopMODE(size_t size_pop, Problem *pro) :
		PopDE<TInd>(size_pop, pro),
		m_mr(0.),
		m_meta(20) {}

	template<typename TInd>
	void PopMODE<TInd>::crossMutate(const std::vector<size_t>& index, TInd& child, Problem *pro, Random *rnd) {
		crossover(this->at(index[0]), this->at(index[1]), this->at(index[2]), child, pro, rnd);
		mutate(child, pro, rnd);
	}

	template<typename TInd>
	void PopMODE<TInd>::crossover(const TInd &parent1, const TInd &parent2, const TInd &parent3, TInd &child, Problem *pro, Random *rnd) {
		int numDim = CAST_CONOP(pro)->numberVariables();
		int idx_rnd = rnd->uniform.nextNonStd(0, numDim);
		auto boundary = CAST_CONOP(pro)->domain();
		Real rate = 0.5;
		child = parent1;
		for (int n = 0; n < numDim; n++) {
			/*Selected Two Parents*/
			child.variable()[n] = parent1.variable()[n] + rate * (parent3.variable()[n] - parent2.variable()[n]);
			if (child.variable()[n] < boundary[n].limit.first) {
				Real r = rnd->uniform.next();
				child.variable()[n] = boundary[n].limit.first + r * (parent1.variable()[n] - boundary[n].limit.first);
			}
			if (child.variable()[n] > boundary[n].limit.second) {
				Real r = rnd->uniform.next();
				child.variable()[n] = boundary[n].limit.second - r * (boundary[n].limit.second - parent1.variable()[n]);
			}
		}

	}

	template<typename TInd>
	void PopMODE<TInd>::mutate(TInd &child, Problem *pro, Random *rnd) {
		int numDim = CAST_CONOP(pro)->numberVariables();
		auto boundary = CAST_CONOP(pro)->domain();
		Real r, delta1, delta2, mut_pow, deltaq;
		Real y, yl, yu, val, xy;
		for (int j = 0; j < numDim; j++) {
			if (rnd->uniform.next() <= m_mr) {
				y = child.variable()[j];
				yl = boundary[j].limit.first;
				yu = boundary[j].limit.second;
				delta1 = (y - yl) / (yu - yl);
				delta2 = (yu - y) / (yu - yl);
				r = rnd->uniform.next();
				mut_pow = 1.0 / (m_meta + 1.0);
				if (r <= 0.5) {
					xy = 1.0 - delta1;
					val = 2.0 * r + (1.0 - 2.0 * r) * (pow(xy, (m_meta + 1.0)));
					deltaq = pow(val, mut_pow) - 1.0;
				}
				else {
					xy = 1.0 - delta2;
					val = 2.0 * (1.0 - r) + 2.0 * (r - 0.5) * (pow(xy, (m_meta + 1.0)));
					deltaq = 1.0 - (pow(val, mut_pow));
				}
				y = y + deltaq * (yu - yl);
				if (y < yl)
					y = yl;
				if (y > yu)
					y = yu;
				child.variable()[j] = y;
			}
		}
	}

	template<typename TInd>
	size_t PopMODE<TInd>::tournamentSelection(Problem *pro, Random *rnd, size_t tournament_size) {
		if (m_rand_seq.size() != this->m_individuals.size()) {
			m_rand_seq.resize(this->m_individuals.size());
			std::iota(m_rand_seq.begin(), m_rand_seq.end(), 0);
		}
		rnd->uniform.shuffle(m_rand_seq.begin(), m_rand_seq.end());
		size_t idx_best = m_rand_seq[0];
		for (size_t k = 1; k < tournament_size; ++k)
			if (this->m_individuals[m_rand_seq[k]]->dominate(*this->m_individuals[idx_best], pro))
				idx_best = m_rand_seq[k];
		return idx_best;
	}

}

#endif //MOEA_DE_H