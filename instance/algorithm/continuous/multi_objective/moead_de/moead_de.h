/********* Begin Register Information **********
{
	"name": "MOEA/D-DE",
	"identifier": "MOEAD_DE",
	"problem tags": [ "MOP", "ConOP" ]
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
// Last modified: 24 Aug 2019 by Xiaofang Wu (email:wuxiaofang@cug.edu.cn)


#ifndef OFEC_MOEAD_DE_H
#define OFEC_MOEAD_DE_H

#include "../../../template/selection/multi_objective/moead.h"
#include "../../../../../core/algorithm/algorithm.h"
#include "../moea_de_pop.h"
#include "../../../../../utility/linear_algebra/vector.h"
#include "../../../../record/multi_objective/rcr_vec_real_moea.h"

namespace ofec {
	class PopMOEAD_DE :public PopMODE<>, MOEAD<IndDE> {
	public:
		PopMOEAD_DE(size_t size_pop, Problem *pro);
		void initialize_(Problem *pro, Random *rnd);
		std::vector<std::pair<Real, Real>> getPopRange() { return m_pop_range; }
		void updatePopRange();
		std::vector<Vector> getWeigh() { return MOEAD::getWeigh(); }
		int evolve(Problem *pro, Algorithm *alg, Random *rnd);
		Population<IndDE>& getOffPop() { return m_off_pop; }

		std::vector<std::vector <std::shared_ptr<Solution<>>>>& getInteractiveSols() { return m_interactive_sol_pair; }
	private:
		std::vector<std::pair<Real, Real>> m_pop_range;
		Population<IndDE> m_off_pop; // offsprings
		//记录交互解
		std::vector<std::vector<std::shared_ptr<Solution<>>>> m_interactive_sol_pair;//2个或3个
	};

	class MOEAD_DE :public Algorithm {
	public:
		void initialize_() override;
		void record() override;
		void initPop();
		std::vector<Vector> getWeigh();
		std::vector<std::pair<Real, Real>> getPopRange() { return m_pop->getPopRange(); }
		void updatePopRange();

		void recordMetrics(Problem *pro, Algorithm *alg);
		std::vector<Real> getIGD() { return m_IGD; }
		void updateHistorySols(Population<IndDE>& pop);
		std::vector <std::shared_ptr<Solution<>>>& getHisSols() { return m_historical_sols; }
		void updateHistoryFrontSols(Population<IndDE>& pop, Problem *pro);
		std::vector <std::shared_ptr<Solution<>>>& getHisFrontSols() { return m_historical_front_sols; }

		


#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;
		std::unique_ptr<PopMOEAD_DE> m_pop;
		size_t m_pop_size;
		std::vector<Real> m_IGD;
		std::vector <std::shared_ptr<Solution<>>> m_historical_sols;//store historical sols
		std::vector <std::shared_ptr<Solution<>>> m_historical_front_sols;//store historical front sols

		

	};

}

#endif
