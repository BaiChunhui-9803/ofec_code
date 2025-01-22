/********* Begin Register Information **********
{
	"name": "NBN_DSA",
	"identifier": "NBN_DynamicSampleAlgorithm",
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


#ifndef NBN_DYNAMIC_SAMPLE_ALGORITHM_H
#define NBN_DYNAMIC_SAMPLE_ALGORITHM_H

#include "../../../core/algorithm/algorithm.h"
#include "../../../utility/nbn_visualization/nbn_grid_division.h"
#include "../../../utility/nbn_visualization/nbn_random_division.h"
#include "../../../utility/nbn_visualization/nbn_visualization_data.h"
#include "../../../utility/nbn_visualization/tree_graph.h"
#include "../../../utility/nbn_visualization/nbn_grid_tree_division.h"
#include "../../../utility/nbn_visualization/nbn_edge_division.h"
#include <memory>




#ifdef  OFEC_DEMO
//#include <custom/buffer/algorithm/combination/point_distribution/buffer_point_distribution.h>
#endif //  OFEC_DEMO


namespace ofec {
	
#define GET_NBN_DYN_SAMPLE(alg) dynamic_cast<NBN_DynamicSampleAlgorithm&>(alg)


	
	class NBN_DynamicSampleAlgorithm : public Algorithm {
	public:

		//const NBN_RandomDivision& nbnNetwork()const {
		//	return m_nbn_network;
		//}
		virtual void record()override {}

		void sample1();
		void samppleCEC2015FO5();
		void sampleTwoObj();

		//void initSolsRand();
		//void initSols2Dgraph();
		void addOpts();


		void outputNBNdata();

		void outputNBNdata(const std::string& path_dir);

		const NBN_VisualizationData& getNBNdata() {
			return m_nbn_data;
		}
		void setMaxSampleNum(int sampleNum) {
			m_maxSample = sampleNum;
		}

		void setFileDir(const std::string& fileDir) {
			m_fileDir = fileDir;
		}
		
	protected:
		virtual void initialize_()override;
		virtual void run_()override;

		void updateBuffer();


		bool isFinish() {
			if (m_nbn_network.m_division->size() < m_minSample) return false;
			else {
				if (m_nbn_network.m_division->size() >= m_maxSample) return true;
				auto curTime = std::chrono::system_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(curTime - m_startTime);
				double duration_seconds = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
				if (duration_seconds >= m_finish_seconds)return true;
				else return false;
			}
		}
	protected:
		///NBN_RandomDivision m_nbn_network;
		std::function<void(SolutionBase& Sol, Problem *pro)> m_eval_fun;
		std::chrono::system_clock::time_point m_startTime;
		double m_duration_seconds = 5;
		double m_finish_seconds = 10*60 ;
//		double m_finish_seconds = 60 ;
		//int m_maxSample = 2e3 ;
		int m_maxSample = 3e6 + 1e3;
		int m_minSample = 3e3;
		std::string m_fileDir;


		std::shared_ptr<SolutionBase> m_cur_sol;
		int m_continousSampleSize;
		int m_initSampleSize = 2e6+1e3;
		const int mc_initSampleSize = 2e3;

		bool  m_filter_flag = false;
		bool  m_opt_flag = false;

		ofec::NBN_DivisionData m_nbn_network;

		NBN_VisualizationData m_nbn_data;
	};

}

#endif 