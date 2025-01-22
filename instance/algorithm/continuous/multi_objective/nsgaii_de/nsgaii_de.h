/********* Begin Register Information **********
{
	"name": "NSGAII-DE",
	"identifier": "NSGAII_DE",
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
// Last modified: 3 May 2019 by Xiaofang Wu (email:wuxiaofang@cug.edu.cn)

#ifndef OFEC_NSGAII_DE_H
#define OFEC_NSGAII_DE_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../moea_de_pop.h"
#include "../../../../record/multi_objective/rcr_vec_real_moea.h"

namespace ofec {
	class PopNSGAII_DE : public PopMODE<>{
	public:
		PopNSGAII_DE(size_t size_pop, Problem *pro);
		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
		Population<IndDE>& getCombinedPop() { return m_pop_combined; }
		Population<IndDE>& getPop() { return *this; }
	protected:
		Population<IndDE> m_pop_combined; // combination of parent and children
	};

	class NSGAII_DE : public Algorithm {
	public:
		void initialize_() override;
		void record() override;
		void initPop();
		void recordMetrics(Problem *pro, Algorithm *alg);

		std::vector<Real> getIGD() { return m_IGD; }
		void updateHistorySols(Population<IndDE>& pop);
		void updateHistoryFrontSols(Population<IndDE>& pop, Problem *pro);
		std::vector <std::shared_ptr<Solution<>>>& getHisSols() { return m_historical_sols; }
		std::vector <std::shared_ptr<Solution<>>>& getHisFrontSols() { return m_historical_front_sols; }

#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;
	protected:
		std::unique_ptr<PopNSGAII_DE> m_pop;
		size_t m_pop_size;
		Real m_mr, m_meta;
		std::vector<Real> m_IGD;
		std::vector <std::shared_ptr<Solution<>>> m_historical_sols;//store historical sols
		std::vector <std::shared_ptr<Solution<>>> m_historical_front_sols;//store historical front sols
	};
}

#endif // !OFEC_NSGAII_DE_H
