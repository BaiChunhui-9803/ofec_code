/********* Begin Register Information **********
{
	"name": "NSGAII-SBX",
	"identifier": "NSGAII_SBX",
	"problem tags": [ "OOMOP","MOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia & Junchen Wang
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 7 Jan 2015
// Last modified: 15 Aug 2019 by Junchen Wang (email:wangjunchen.chris@gmail.com)
// Last modified: 23 Aug 2019 by Xiaofang Wu (email:email:wuxiaofang@cug.edu.cn)

#ifndef OFEC_NSGAII_SBX_H
#define OFEC_NSGAII_SBX_H

#include "../../../template/classic/ga/sbx_pop.h"
//#include "../../../template/selection/multi_objective/nsgaii.h"
#include "../../../../../core/algorithm/algorithm.h"
#include "../../../../record/multi_objective/rcr_vec_real_moea.h"

namespace ofec {
	class PopNSGAII_SBX : public PopSBX<>{
	public:
		PopNSGAII_SBX(size_t size_pop, Problem *pro);
		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
		Population<Solution<>>& getCombinedPop() { return m_pop_combined; }
		Population<Solution<>>& getPop() { return *this; }
		Population<Solution<>>& getOffspring() { return m_offspring; }
	protected:
		Population<Solution<>> m_pop_combined;  // combination of parent and children
		Population<Solution<>> m_offspring;  // 
	};

	class NSGAII_SBX : public Algorithm {
	protected:
		std::unique_ptr<PopNSGAII_SBX> m_pop;
		size_t m_pop_size;
		Real m_cr, m_mr, m_ceta, m_meta;
		std::vector<Real> m_IGD;
		std::vector <std::shared_ptr<Solution<>>> m_historical_sols;//store historical sols
		std::vector <std::shared_ptr<Solution<>>> m_historical_front_sols;//store historical front sols

		void initialize_() override;
		void run_() override;
		void initPop();
		
		std::vector<Real> getIGD() { return m_IGD; }
		void updateHistorySols(Population<Solution<>>& pop);
		std::vector <std::shared_ptr<Solution<>>>& getHisSols() { return m_historical_sols; }
		void updateHistoryFrontSols(Population<Solution<>>& pop, Problem *pro);
		std::vector <std::shared_ptr<Solution<>>>& getHisFrontSols() { return m_historical_front_sols; }

#ifdef OFEC_DEMO
		void updateBuffer();
#endif

	public:
		void record() override;
		void recordMetrics(Problem *pro, Algorithm *alg);
	};
}

#endif // !OFEC_NSGAII_SBX_H

