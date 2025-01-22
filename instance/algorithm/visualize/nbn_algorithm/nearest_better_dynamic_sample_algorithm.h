/********* Begin Register Information **********
{
	"name": "VISUAL_NBN_SAMPLE",
	"identifier": "NearestBetterDynamicSampleAlgorithm",
	"problem tags": [ 
		"SOP", "MOP", "DOP", "DMOP", "MMOP", "GOP", "ROOT", 
		"ConOP", "ComOP", "TSP", "COP", "VRP", "TTP", "JSP", 
		"KOP", "SAT", "OneMax", "QAP", "MKP", "EOP", "LSOP",
		"CSIWDN", "DVRP", "SP", "APP", "NoisyOP", "PMP"
	]
}
*********** End Register Information **********/


/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@163.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/


#ifndef NEAREST_BETTER_DYNAMIC_SAMPLE_ALGORITHM_H
#define NEAREST_BETTER_DYNAMIC_SAMPLE_ALGORITHM_H

#include "../../../core/algorithm/algorithm.h"

#include "nearest_better_network_sample.h"
#include "sample2d_graph_algorithm.h"

#include <memory>



#ifdef  OFEC_DEMO
//#include <custom/buffer/algorithm/combination/point_distribution/buffer_point_distribution.h>
#endif //  OFEC_DEMO


namespace ofec {
#define GET_NBNSample(alg) dynamic_cast<NearestBetterDynamicSampleAlgorithm&>(alg)
	class NearestBetterDynamicSampleAlgorithm : public Algorithm {
	public:

		const NearestBetterNetworkSample& nbnNetwork()const {
			return m_nbn_network;
		}
		virtual void record()override {}

		void sample1();
		void samppleCEC2015FO5();
		void sampleTwoObj();

		void initSolsRand();
		void initSols2Dgraph();
		void addOpts();
	protected:
		virtual void initialize_()override;
		virtual void run_()override;

		void updateBuffer();
	protected:
		NearestBetterNetworkSample m_nbn_network;
		std::function<void(SolutionBase& Sol, Problem *pro)> m_eval_fun;
		double m_duration_seconds = 5;

		std::shared_ptr<SolutionBase> m_cur_sol;
		Sample2D_Graph_Algorithm m_sample_alg;
		int m_continousSampleSize;
		int m_initSampleSize = 2e3;
		const int mc_initSampleSize = 2e3;

		bool  m_filter_flag = false;
		bool  m_opt_flag = false;
	};

}

#endif 