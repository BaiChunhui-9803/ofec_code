/********* Begin Register Information **********
{
	"name": "VISUAL_DSA",
	"identifier": "DynamicSampleAlgorithm",
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


#ifndef DYNAMIC_SAMPLE_ALGORITHM_H
#define DYNAMIC_SAMPLE_ALGORITHM_H

#include "../../../core/algorithm/algorithm.h"
#include "sample_nearest_better_network.h"
#include "sample2d_graph_algorithm.h"
#include "../../../core/problem/solution.h"



#ifdef  OFEC_DEMO
//#include <custom/buffer/algorithm/combination/point_distribution/buffer_point_distribution.h>
#endif //  OFEC_DEMO


namespace ofec {
#define GET_NBN(alg) dynamic_cast<DynamicSampleAlgorithm&>(alg)


	class DynamicSampleAlgorithm : public Algorithm {
	public:
		using SolutionType = typename Solution<> ;
	public:
		virtual void record()override;
		const SampleNearestBetterNetwork& nbnNetwork()const {
			return m_nbn_network;
		}
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
		SampleNearestBetterNetwork m_nbn_network;
		std::function<void(Solution<>& sol, Problem *pro)> m_eval_fun;
		double m_duration_seconds = 5;

		std::shared_ptr<SolutionType> m_cur_sol;
		Sample2D_Graph_Algorithm m_sample_alg;
		int m_continousSampleSize;
		int m_initSampleSize = 2e3;
		const int mc_initSampleSize = 2e3;

		bool  m_filter_flag = false;
		bool  m_opt_flag = false;
	};

}

#endif 