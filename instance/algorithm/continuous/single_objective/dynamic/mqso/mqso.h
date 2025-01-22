/********* Begin Register Information **********
{
	"name": "mQSO",
	"identifier": "mQSO",
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
// Created by PengMai at 5 July 2021
// Updated by PengMai at 14 August 2021 

/* ---------------------------------------------------------------------------------------
mQSO & mCPSO: Blackwell T, Branke J.
Multiswarms, exclusion, and anti-convergence in dynamic environments[J].
IEEE Transactions on Evolutionary Computation, 2006, 10(4):459-472.
-----------------------------------------------------------------------------------------*/
#ifndef OFEC_MQSO_H
#define OFEC_MQSO_H

#include "mqso_subswarm.h"
#include "../metrics_dynamic.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../problem/continuous/single_objective/dynamic/uncertianty_continuous.h"
#include "../../../../../record/dynamic/rcr_vec_real_dynamic.h"

namespace ofec {
	class mQSO : public MetricsDynamicConOEA {
	protected:
		void initialize_() override;
		void initializeSolution(int num);
		bool checkConvergenceAll();
		void exclude();
		void evolve();
		void run_() override;
		int updateBestSubPopFlag();
		void updateSubPopBestIdx();
		int findWorstPop();
		void autiConvergence();
		void removeRedundancePop();
		void measureMultiPop(bool);
	public:
		void record() override;
#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		MultiPopulation<MQSOSwarm> m_sub_pop;
		int m_M;        // the number of populations
		double m_Rconv; // threshold radius of conveged swarms
		double m_Rexcl; // radius of exlusion radius 
		bool m_exclusion;	// 1 for with exclusion and 0 no exclustion	
		int mc_init_indiSize;

		bool m_convergedAll;
		int m_bestsubPopIdx;

		long long int count_converge = 0;
		long long int count_exclusion = 0;
	};

	using mCPSO = mQSO;		//for mCPSO -> charged particle
}

#endif //!OFEC_MQSO_H