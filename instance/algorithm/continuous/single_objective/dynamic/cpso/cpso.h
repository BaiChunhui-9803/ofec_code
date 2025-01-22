/********* Begin Register Information **********
{
	"name": "CPSO",
	"identifier": "CPSO",
	"problem tags": [ "ConOP", "SOP", "GOP", "MMOP", "DOP" ]
}
*********** End Register Information **********/

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
// Created by Changhe Li
// Updated by Mai Peng at 26 September 2021 

/* ---------------------------------------------------------------------------------------
S. Yang and C. Li,
A clustering particle swarm optimizer for locating and tracking multiple optima in dynamic environments,
IEEE Trans. Evol. Comput., vol. 14, no. 6, pp. 959-974, 2010.
-----------------------------------------------------------------------------------------*/

#ifndef OFEC_CPSO_H
#define OFEC_CPSO_H
#include "cpso_subswarm.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../problem/continuous/single_objective/dynamic/uncertianty_continuous.h"
#include "../../../../../record/dynamic/rcr_vec_real_dynamic.h"
#include "../../../../../../utility/clustering/hslh.h"
#include "../metrics_dynamic.h"
#ifdef OFEC_DEMO
#include "../dynamic_pso.h"
#endif

namespace ofec {
	class CPSO : public MetricsDynamicConOEA
#ifdef OFEC_DEMO
		, public DynamicPSO
#endif
	{
	protected:
		void initialize_() override;
		int initializeOriginalPop(int num);
		int createSubswarms(int num);
		int checkOverlapping();
		int checkConverging();

		void run_() override;

		void informChange(int rf);
		void measureMultiPop(bool flag);

	public:
		void record() override;
#ifdef OFEC_DEMO
		void updateBuffer();
		std::vector<bool> getPopHiberState() const override;
#endif
	protected:
		ofec::MultiPopulation<CPSOSwarm> m_sub_pop;					//slst
		int m_init_popsize;
		int m_min_subpopSize = 1;
		std::vector<std::unique_ptr<templateParticle>> m_all_indis;			//initial cradle C
		std::vector<templateParticle> m_clst_indis;		//clst
		int m_max_subpopSize;
		Real m_overlapDegree;
		Real m_rconv;
	};
}
#endif //!OFEC_CPSO_H
