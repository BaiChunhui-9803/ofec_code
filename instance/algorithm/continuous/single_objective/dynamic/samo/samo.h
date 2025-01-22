/********* Begin Register Information **********
{
	"name": "SAMO",
	"identifier": "SAMO",
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
// Created by PengMai at 20 July 2021
// Updated by PengMai at 16 August 2021 

/* ---------------------------------------------------------------------------------------
SAMO: Blackwell, Tim & Branke, Juergen. (2004).
Particle Swarm Optimization in Dynamic Environments. 3005. 489-500.
-----------------------------------------------------------------------------------------*/

#ifndef OFEC_SAMO_H
#define OFEC_SAMO_H

#include "../mqso/mqso_subswarm.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include "../metrics_dynamic.h"
#include "../../../../../record/dynamic/rcr_vec_real_dynamic.h"

namespace ofec {
	class SAMO : public MetricsDynamicConOEA {
	public:
		void initialize_();
		void initializeSolution(int num);
		int findWorstPop();
		void checkSwarmConverge();
		void checkSwarmConverge2();
		void checkSwarmState();
		int reInitializePop(int idx);
		void removeRedundancePop();
		void evolve();
		void exclude();
		void run_() override;
		void record() override; 
		void measureMultiPop(bool);
		void updateRexcludeRaduis();
#ifdef OFEC_DEMO
		void updateBuffer();
		int getMfree() const { return m_Mfree; }
#endif
	protected:
		MultiPopulation<MQSOSwarm> m_sub_pop;
		//		int m_M;        // remove from mQSO & mCPSO
		double m_Rconv; // threshold radius of conveged swarms
		double m_Rexcl; // radius of exlusion radius 
		bool m_exclusion;	// 1 for exclusion and 0 no exclustion

		int m_M;
		int m_Mfree;
		int m_Necxess;
		long long int count_converge = 0;
		long long int count_exclusion = 0;
	};
}
#endif