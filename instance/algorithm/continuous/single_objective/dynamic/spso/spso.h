/********* Begin Register Information **********
{
	"name": "SPSO",
	"identifier": "SPSO",
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
// Created by PengMai at 28 Jun 2021
// Updated by PengMai at 13 September 2021 

/* ---------------------------------------------------------------------------------------
SPSO: [1] Parrott, D.; Xiaodong Li; , "Locating and tracking multiple dynamic optima
by a particle swarm model using speciation," Evolutionary Computation, IEEE Transactions
on , vol.10, no.4, pp.440-458, Aug. 2006.
-----------------------------------------------------------------------------------------*/
//? the result isn't good as the paper

#ifndef OFEC_SPSO_H
#define OFEC_SPSO_H
#include "spso_subswarm.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../problem/continuous/single_objective/dynamic/uncertianty_continuous.h"
#include "../../../../../record/dynamic/rcr_vec_real_dynamic.h"
#include "../metrics_dynamic.h"
#ifdef OFEC_DEMO
#include "../dynamic_pso.h"
#endif

namespace ofec {
	class SPSO : public MetricsDynamicConOEA
#ifdef OFEC_DEMO
		, public DynamicPSO
#endif
	{
	protected:

		void initialize_() override;
		void initializeOriginalPop(int num);
		void sortParticles();
		void informChange(int rf);
		void measureMultiPop(bool flag);
		void measureAllParticles(bool flag);
		void getAllParticles();
		void createSpecies();
		void updateReplaceParticles();
		bool removeConvergence();
		void checkOvercrowd();

		void run_() override;

	public:
		void record() override;
#ifdef OFEC_DEMO
		void updateBuffer();
		std::vector<bool> SPSO::getPopHiberState() const override;
#endif
	protected:
		ofec::MultiPopulation<SPSOSwarm> m_sub_pop;
		int m_all_indisize;
		int m_num_seeds = 0;
		std::vector<int> m_seedList;
		int m_pMax = 10;
		std::vector<SPSOParticle> m_all_indis;
		std::vector<SPSOParticle> m_remain_indis;
		std::vector<SPSOParticle> m_replace_indis;
		std::vector<SPSOParticle> m_seed;
		std::vector<std::vector<SPSOParticle>> m_species;

		Real m_rs;
		bool first_ajustment = true;

	};
}
#endif //!OFEC_SPSO_H
