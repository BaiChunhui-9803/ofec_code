/********* Begin Register Information **********
{
	"name": "HmSO",
	"identifier": "HmSO",
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
// Updated by PengMai at 5 July 2021 

/* ---------------------------------------------------------------------------------------
HmSO: [1] Kamosi, M. ,  A. B. Hashemi , and  M. R. Meybodi .
A hibernating multi-swarm optimization Algorithm for dynamic environments.
2010 Second World Congress on Nature and Biologically Inspired Computing (NaBIC) IEEE, 2011.
-----------------------------------------------------------------------------------------*/

#ifndef OFEC_HMSO_H
#define OFEC_HMSO_H
#include "hmso_subswarm.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../problem/continuous/single_objective/dynamic/uncertianty_continuous.h"
#include "../../../../../record/dynamic/rcr_vec_real_dynamic.h"
#include "../metrics_dynamic.h"
#ifdef OFEC_DEMO
#include "../dynamic_pso.h"
#endif

namespace ofec {
	class HmSO : public MetricsDynamicConOEA
#ifdef OFEC_DEMO
		, public DynamicPSO
#endif
	{
	protected:
		void initialize_() override;
		void initializeParent(int num);
		void initializeChild(HmSOSwarm& t, int num);
		void checkOverlapping();
		void removeOverlapping();
		void createChildSwarm();
		void updateChildSwarm();
		int updateBestSubPopIdx();
		void updateChildAttractor();
		void wakeupHiberPop();
		void relaunchChildSwarm(HmSOSwarm&);
		void informChange(int rf);
		bool checkParentConvergence();
		void preventParentConvergence();
		void measureMultiPop(bool flag);
		void run_() override;

	public:
		void record() override;
#ifdef OFEC_DEMO
		void updateBuffer();
		std::vector<bool> getPopHiberState() const override;
#endif
	protected:
		ofec::MultiPopulation<HmSOSwarm> m_sub_pop;
		Real m_rpc;	//radius between parent & child
		Real m_rconv; // threshold radius of conveged swarms
		Real m_rexcl; // radius of exlusion radius
		Real m_rs; // radius of random search
		Real m_margin; // margin for hibernation
		std::unique_ptr<Solution<>> m_preiter;
		std::unique_ptr<Solution<>> m_curiter;
	};


}
#endif //!OFEC_HMSO_H
