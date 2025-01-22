/********* Begin Register Information **********
{
	"name": "MOEA/D-SBX",
	"identifier": "MOEAD_SBX",
	"problem tags": [ "ConOP", "MOP" ]
}
*********** End Register Information **********/

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

*************************************************************************/
// Created: 13 Jan 2015
// Last modified: 25 Aug 2019 by Xiaofang Wu (email:wuxiaofang@cug.edu.cn)


#ifndef OFEC_MOEAD_SBX_H
#define OFEC_MOEAD_SBX_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../../../template/selection/multi_objective/moead.h"
#include "../../../template/classic/ga/sbx_pop.h"
#include "../../../../record/multi_objective/rcr_vec_real_moea.h"

namespace ofec {
	class PopMOEAD_SBX :public PopSBX<>, MOEAD<Solution<>> {
	public:
		PopMOEAD_SBX(size_t size_pop, Problem *pro);
		void initialize_(Problem *pro, Random *rnd);
		std::vector<std::pair<Real, Real>> getPopRange() { return m_pop_range; }
		void updatePopRange();
		std::vector<Vector> getWeigh() { return MOEAD::getWeigh(); }
		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
		Population<Solution<>>& getOffPop() { return m_off_pop; }
	private:
		std::vector<std::pair<Real, Real>> m_pop_range;
		Population<Solution<>> m_off_pop; // offsprings
	};


	class MOEAD_SBX :public Algorithm {
	public:
		void initialize_() override;
		void record() override;
		void initPop();
		std::vector<Vector> getWeigh();
		std::vector<std::pair<Real, Real>> getPopRange() { return m_pop->getPopRange(); }
		void updatePopRange();

		void recordMetrics(Problem *pro, Algorithm *alg);
		std::vector<Real> getIGD() { return m_IGD; }
		void updateHistorySols(Population<Solution<>>& pop);
		std::vector <std::shared_ptr<Solution<>>>& getHisSols() { return m_historical_sols; }
		void updateHistoryFrontSols(Population<Solution<>>& pop, Problem *pro);
		std::vector <std::shared_ptr<Solution<>>>& getHisFrontSols() { return m_historical_front_sols; }

#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;
		std::unique_ptr<PopMOEAD_SBX> m_pop;
		size_t m_pop_size;
		Real m_cr, m_mr, m_ceta, m_meta;

		std::vector<Real> m_IGD;
		std::vector <std::shared_ptr<Solution<>>> m_historical_sols;//store historical sols
		std::vector <std::shared_ptr<Solution<>>> m_historical_front_sols;//store historical front sols

	};

}

#endif //!OFEC_MOEAD_SBX_H
