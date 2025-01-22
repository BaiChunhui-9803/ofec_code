/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*-------------------------------------------------------------------------------
*
*  Modified by Diao Yiya on 2024 09 16
*
*-------------------------------------------------------------------------------
*************************************************************************/


#ifndef OFEC_SAMPLING_ALGORITHM_MULTIPOP_MULTITHREAD_H
#define OFEC_SAMPLING_ALGORITHM_MULTIPOP_MULTITHREAD_H


#include <vector>

#include "sampling_algorithm_multipop.h"
#include "../../../../core/environment/environment.h"
#include "../../../../utility/general_multithread/general_multithread.h"
#include "../../../../run/interface.h"





namespace ofec {


	struct SamplingAgorithmMutithread {

		class RunAlgTaskLocalThreadInfo : public GeneralMultiThreadInfo {
		public:
			double m_seed = 0.1;
			int m_algId = 0;
			std::string m_algName;
			ofec::ParameterMap m_algParam;
			//	double m_firstRandom = 0;
		};

		class RunAlgTaskGlobalThreadInfo : public GeneralMultiThreadInfo {
		public:

			double m_pro_seed = 0.5;
			double m_env_seed = 0.5;
			std::string m_proName;
			ofec::ParameterMap m_proParam;

			std::vector<std::vector<std::shared_ptr<ofec::SolutionBase>>> m_sols;
			std::vector<std::vector<std::shared_ptr<SamplingData::SolutionInfo>>> m_solInfos;
		};


		static void runAlgTask(
			std::unique_ptr<GeneralMultiThreadInfo>& curThreadInfo,
			std::unique_ptr<GeneralMultiThreadInfo>& threadLocalInfo,
			std::unique_ptr<GeneralMultiThreadInfo>& globalThreadInfo) {

			auto& localInfo = dynamic_cast<RunAlgTaskLocalThreadInfo&>(*curThreadInfo);
			//	auto& threadInfo = dynamic_cast<ThreadLocalInfo&>(*threadLocalInfo);
			auto& globalInfo = dynamic_cast<RunAlgTaskGlobalThreadInfo&>(*globalThreadInfo);

			std::shared_ptr<Environment> env_ptr(Environment::create());
			auto& env = *env_ptr;

			env.initialize(globalInfo.m_env_seed);


			env.setProblem(Factory<Problem>::produce(globalInfo.m_proName));
			env.problem()->inputParameters().input(globalInfo.m_proParam);
			env.problem()->recordInputParameters();
			env.initializeProblem(globalInfo.m_pro_seed);
			env.setAlgorithm(Factory<Algorithm>::produce(localInfo.m_algName));
			env.algorithm()->inputParameters().input(localInfo.m_algParam);
			env.algorithm()->recordInputParameters();
			env.initializeAlgorithm(localInfo.m_seed);
			env.runAlgorithm();
			
			auto samplingData = CAST_SAMPLE_DATA(env.algorithm());
			globalInfo.m_sols[localInfo.m_algId] = samplingData->getSols();
			globalInfo.m_solInfos[localInfo.m_algId] = samplingData->getSolInfos();
		}

		static void
			runAlgMultiTask(
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
				std::vector<std::shared_ptr<SamplingData::SolutionInfo>>& solInfos,
				double envSeed,
				double proSeed,
				const std::string& proName,
				const  ofec::ParameterMap& pro_param,
				std::string& algName,
				const ofec::ParameterMap& algParam,
				int numRun) {
			ofec::GeneralMultiThread curThread;
			curThread.ms_global_info.reset(new RunAlgTaskGlobalThreadInfo());
			auto& globalInfo(dynamic_cast<RunAlgTaskGlobalThreadInfo&>(*curThread.ms_global_info));
			globalInfo.m_pro_seed = proSeed;
			globalInfo.m_env_seed = envSeed;
			globalInfo.m_proName = proName;
			globalInfo.m_proParam = pro_param;
			globalInfo.m_sols.resize(numRun);
			globalInfo.m_solInfos.resize(numRun);
			//globalInfo.m_firstRand.resize(numRun);
			//globalInfo.m_lastRand.resize(numRun);
			RunAlgTaskLocalThreadInfo localInfo;

			for (int idx(0); idx < numRun; ++idx) {
				localInfo.m_algName = algName;
				localInfo.m_algParam = algParam;
				localInfo.m_algId = idx;
				localInfo.m_seed = double(idx + 1) / double(numRun + 1);
				curThread.ms_info_buffer.emplace_back(new RunAlgTaskLocalThreadInfo(localInfo));
			}

			curThread.ms_fun = std::bind(&SamplingAgorithmMutithread::runAlgTask, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
			// out(totalInfo.m_save_dir + "error.txt");
			curThread.ms_cout_flag = false;
			curThread.setDefaultNumberThread();
			curThread.runNoException();

			for (auto& it : globalInfo.m_sols) {
				for (auto& it2 : it) {
					sols.emplace_back(it2);
				}
			}

			for (int runId(0); runId < numRun; ++runId) {
				for (auto& it2 : globalInfo.m_solInfos[runId]) {
					it2->m_runId = runId;
					solInfos.emplace_back(it2);
				}
			}
		}


	};
}

#endif