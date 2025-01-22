/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*********************************************************************************/
// updated Jan 21, 2019 by Li Zhou

#ifndef OFEC_CC_FRAMEWORK_H
#define OFEC_CC_FRAMEWORK_H

#include <memory>
#include "../../../../core/algorithm/population.h"
#include "../../../../utility/functional.h"

namespace ofec {
	template <typename TPop1, typename TPop2, typename TInd>
	class CC_framework : public Population<TInd>   // just for epanet
	{

	public:
		enum class FillType { random, best };
		std::unique_ptr<TPop1> m_first_pop;
		std::unique_ptr<TPop2> m_second_pop;
		std::unique_ptr<TInd> m_best;

		std::vector<std::vector<float>> m_history_mass;

	public:
		CC_framework() = default;
		CC_framework(size_t size,Environment *env) {
			auto pro = env->problem();
			m_first_pop = std::move(std::unique_ptr<TPop1>(new TPop1(size, env, CAST_CSIWDN(pro)->numSource())));
			m_second_pop = std::move(std::unique_ptr<TPop2>(new TPop2(size, env, CAST_CSIWDN(pro)->numSource())));
			size_t num_obj = CAST_CSIWDN(pro)->numberObjectives();
			size_t num_cons = CAST_CSIWDN(pro)->numberConstraints();
		}
		virtual ~CC_framework() {}

		const TPop1& first_pop() const{
			return *m_first_pop;
		}
		
		const TPop2& second_pop()const {
			return *m_second_pop;
		}


		int evolve(Environment *env, Random *rnd) override { return kNormalEval; }
		int evolve(Environment *env, Random *rnd, const std::pair<int, int>& source_index);
		int evaluate(Environment *env) override;
		
		CC_framework& operator=(CC_framework & cc) {
			*m_first_pop = *cc.m_first_pop;
			*m_second_pop = *cc.m_second_pop;
			if (m_best)
				*m_best = *cc.m_best;
			else
				m_best.reset(new TInd(*cc.m_best));
			return *this;
		}

		CC_framework& operator=(CC_framework && cc) {
			*m_first_pop = std::move(*cc.m_first_pop);
			*m_second_pop = std::move(*cc.m_second_pop);
			if (m_best)
				*m_best = *cc.m_best;
			else
				m_best.reset(cc.m_best.release());
			return *this;
		}

		int fillEach(Random *rnd, Environment *env, FillType type, const std::pair<int, int>& source_index);
		void updateBest(Environment *env);

		void initializeNewPop(Environment *env, Random *rnd, VarCSIWDN& var_best, const std::pair<int, int>& source_index);
		bool isStableMass();
		void recordHistoryMass(Environment *env, const std::pair<int, int>& source_index);
	};

	template<typename TPop1, typename TPop2, typename TInd>
	int CC_framework<TPop1, TPop2, TInd>::evolve(Environment *env, Random *rnd, const std::pair<int, int> &source_index) {

		m_first_pop->evolve(env, rnd, false, source_index);
		
		m_second_pop->updateCR(env,rnd, source_index);
		m_second_pop->updateF(env, rnd, source_index);
		m_second_pop->updateProStrategy(env, source_index);
		m_second_pop->mutate(env, rnd, source_index);
		m_second_pop->recombine(env, rnd, source_index);
		m_second_pop->select(env,  false, source_index);

		m_first_pop->increaseIteration();
		m_second_pop->increaseIteration();
		
		this->increaseIteration();
		return fillEach(rnd, env, CC_framework<TPop1, TPop2, TInd>::FillType::best, source_index);

	}

	template<typename TPop1, typename TPop2, typename TInd>
	int CC_framework<TPop1, TPop2, TInd>::evaluate(Environment *env) {
		int tag = kNormalEval;

		for (size_t i = 0; i < m_first_pop->size(); ++i) {
			tag = m_first_pop->at(i).evaluate(env);
			if (tag != kNormalEval) return tag;
		}
		tag = m_first_pop->best().front()->evaluate(env);
		if (tag != kNormalEval) return tag;

		for (size_t i = 0; i < m_second_pop->size(); ++i) {
			tag = m_second_pop->at(i).evaluate(env);
			if (tag != kNormalEval) return tag;
		}
		tag = m_second_pop->best().front()->evaluate(env);
		if (tag != kNormalEval) return tag;

		tag = m_best->evaluate(env);

		return tag;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	int CC_framework<TPop1, TPop2, TInd>::fillEach(Random *rnd, Environment *env, FillType type, const std::pair<int, int> &source_index) {
		if (type == FillType::random) {
			size_t I = rnd->uniform.nextNonStd<size_t>(0, m_second_pop->size());
			m_first_pop->fillSolution((*m_second_pop)[I].variable(), env, source_index);
			I = rnd->uniform.nextNonStd<size_t>(0, m_second_pop->size());
			m_second_pop->fillSolution((*m_first_pop)[I].variable(), env, source_index);
		}
		else if (type == FillType::best) {
			m_first_pop->updateBest(env);
			m_second_pop->updateBest(env);
			m_first_pop->fillSolution(m_second_pop->best().front()->variable(), env, source_index);
			m_second_pop->fillSolution(m_first_pop->best().front()->variable(), env, source_index);
		}
		m_first_pop->evaluate(env);
		m_second_pop->evaluate(env);

		//updateBest(pro, alg);
		//m_history_best_obj.push_back(m_best->objective()[0]);
		return kNormalEval;
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void CC_framework<TPop1, TPop2, TInd>::updateBest(Environment *env) {
		auto pro = env->problem();
		m_first_pop->updateBest(env);
		m_second_pop->updateBest(env);
		auto &best = dominate(*m_first_pop->best().front(),*m_second_pop->best().front(), pro->optimizeMode()) ?
			*m_first_pop->best().front() : *m_second_pop->best().front();
		if (!m_best)
			m_best.reset(new TInd(best));
		else if (dominate(best, *m_best, pro->optimizeMode()))
			*m_best = best;



		//size_t number_objectives = pro->numberObjectives();
		//size_t num_cons = pro->numberConstraints();
		//size_t dim = CAST_CSIWDN(pro)->numberSource();
		//TInd temp(number_objectives, num_cons, dim);

		//int z = 0;
		//while ((z + 1) < CAST_CSIWDN(pro)->numberSource() && CAST_CSIWDN(pro)->phase() >= ((*(m_first_pop->best().front())).variable().startTime(z + 1) / CAST_CSIWDN(pro)->intervalTimeStep())) {
		//	z++;
		//}

		//temp = TInd(*(m_first_pop->best().front()));
		///*temp.variable().duration(i) = (*(m_first_pop->best().front())).variable().duration(i);
		//temp.variable().startTime(i) = (*(m_first_pop->best().front())).variable().startTime(i);
		//temp.variable().multiplier(i) = (*(m_first_pop->best().front())).variable().multiplier(i);*/

		//temp.evaluate(pro, alg);

		//if (m_best->variable().index(z) == -1) {
		//	*m_best = temp;
		//}
		//else {
		//	m_best->evaluate(pro, alg);
		//	if (temp.objective()[0] < m_best->objective()[0])
		//		*m_best = temp;
		//}
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void CC_framework<TPop1, TPop2, TInd>::initializeNewPop(Environment *env, Random *rnd, VarCSIWDN& var_best, const std::pair<int, int>& source_index) {
		m_first_pop->initialize(env, rnd);
		m_second_pop->initialize(env, rnd);
		int q = source_index.first, z = source_index.second;
		for (size_t i = 0; i < m_first_pop->size(); i++) {
			for (size_t j = 0; j < z; j++) {
				(*m_first_pop)[i].variable().getEpa(j) = var_best.getEpa(j);
				(*m_second_pop)[i].variable().getEpa(j) = var_best.getEpa(j);
				(*m_first_pop)[i].mpu().variable().getEpa(j) = var_best.getEpa(j);
				(*m_first_pop)[i].mpv().variable().getEpa(j) = var_best.getEpa(j);
				(*m_second_pop)[i].mpu().variable().getEpa(j) = var_best.getEpa(j);
				(*m_second_pop)[i].mpv().variable().getEpa(j) = var_best.getEpa(j);
			}
		}
		m_first_pop->evaluate(env);
		m_first_pop->updateBest(env);
		m_second_pop->evaluate(env);
		m_second_pop->updateBest(env);
		updateBest(env);
	}

	template <typename TPop1, typename TPop2, typename TInd>
	void CC_framework<TPop1, TPop2, TInd>::recordHistoryMass(Environment *env, const std::pair<int, int>& source_index) {
		int z = source_index.second;
		m_history_mass.push_back(m_best->variable().multiplier(z));
	}

	template <typename TPop1, typename TPop2, typename TInd>
	bool CC_framework<TPop1, TPop2, TInd>::isStableMass() {
		size_t standard = 3;
		if (m_history_mass.size() < standard) return false;
		int i = m_history_mass.size() - 1;
		int count = 0;
		while (m_history_mass[i] == m_history_mass[i - 1]) {
			++count;
			--i;
			if (i == 0) break;
		}
		return count >= standard - 1;
	}
}

#endif