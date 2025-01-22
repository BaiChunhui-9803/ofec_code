/********* Begin Register Information **********
{
	"name": "DNSGAII-A",
	"identifier": "DNSGAII_A",
	"problem tags": [ "DMOP", "MOP", "ConOP", "DOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Qingshan Tan
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

* see https://github.com/Changhe160/OFEC for more information
*************************************************************************/

/************************************************************************
Deb K, Karthik S.
Dynamic multi-objective optimization and decision-making using modified NSGA-II:
a case study on hydro-thermal power scheduling[C]
International conference on evolutionary multi-criterion optimization. Springer, Berlin, Heidelberg, 2007: 803-817.
************************************************************************/

// Created: 19 Octorber 2019 by Qingshan Tan

#ifndef OFEC_DNSGAII_A_H
#define OFEC_DNSGAII_A_H

//#include "../../../../../template/selection/multi_objective/nsgaii.h"
#include "../../../template/classic/ga/sbx_pop.h"
#include "../../../../../core/algorithm/algorithm.h"
#include "../../../../../core/algorithm/population.h"
#include "../detector.h"

namespace ofec {
	class DNSGAII_A_pop : public PopSBX<>{
	public:
		DNSGAII_A_pop(size_t size_pop, Problem* pro, size_t size_var, size_t size_obj, const std::vector<OptimizeMode>& opt_mode);
		void initialize_(Problem* pro, Random* rnd);
		int evolve(Problem* pro, Algorithm* alg, Random* rnd) override;
		//std::vector<size_t> get_rand_seq() { return m_rand_seq; }
		
		bool problemChanged(Real r, Problem* pro, Algorithm* alg, Random* rnd);
		bool populationConverged();
		void responseChange(Problem* pro, Algorithm* alg, Random* rnd);
		void updatePop(Problem* pro, Random* rnd);
		void addDiversity(Problem* pro, Algorithm* alg, Random* rnd);

		Real getDetectRate() { return m_detect_rate; }
		Real getReplaceRate() { return m_replace_rate; }
		void setDetectRate(Real r) { m_detect_rate = r; }
		void setReplaceRate(Real r) { m_replace_rate = r; }
	protected:
		std::vector<Solution<>> m_offspring;  // 2 size of population
		Real m_detect_rate = 0.1;//the ratio of population to re-evaluate
		Real m_replace_rate = 0.2;//the ratio of population to be replaced
	};

	class DNSGAII_A : public Algorithm {
	public:
		void initialize_() override;
		void record() override;
		void initPop();

#ifdef OFEC_DEMO
		void updateBuffer() {}
#endif
	protected:
		void run_() override;
	protected:
		std::unique_ptr<DNSGAII_A_pop> m_pop;
		size_t m_pop_size;
		Real m_cr, m_mr, m_ceta, m_meta;
		
	};
}

#endif // !OFEC_DNSGAII_A_H
