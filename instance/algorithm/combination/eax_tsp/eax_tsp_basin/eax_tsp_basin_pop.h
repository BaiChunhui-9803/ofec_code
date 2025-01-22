/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*-------------------------------------------------------------------------------
*
*  Created by Diao Yiya on 2023 11 13
*
*-------------------------------------------------------------------------------
*

*
*************************************************************************/

#ifndef OFEC_EAX_TSP_BASIN_POP_H
#define OFEC_EAX_TSP_BASIN_POP_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/problem/solution.h"
#include "../eax_tsp_multithread/environment.h"
#include "hnsw_basin.h"
#include <mutex>


namespace ofec {



	class PopEaxTspBasin : public Population<Solution<VariableVector<int>>>,
		public eax_tsp_mt::TEnvironment {
		using SolutionType = eax_tsp_mt::TIndi;
		

		struct ThreadInfo {
			std::shared_ptr<Random> m_random;
			std::shared_ptr<eax_tsp_mt::TCross> tCross;         /* Eede assembly crossover */
			std::shared_ptr<eax_tsp_mt::TKopt> tKopt;           /* Local search with the 2-opt neighborhood */
		};

	protected:

		int m_curBasinId = -1;
		int m_pos = 0;
		int m_curIter = 0;
		int m_maxLoopInGeneration = 1e3;
		int m_numThread = 0;
		std::vector<ThreadInfo> m_threadInfos;
		//std::mutex m_mutex;

		bool m_terminatedFlag = false;
		
		std::vector<std::shared_ptr<SolutionType>> m_history_data;
		std::vector<SolutionType*> m_curPop;
		SolutionType* m_best = nullptr;

		void initThreadInfo(ThreadInfo& cur, Random* rnd);
		
		
		void updateFitness(SolutionType* indi) {
			indi->setFitness(m_pos * indi->fEvaluationValue);
		}
		bool judgeInsideBasin(const HNSWbasin& model, SolutionBase& sol, Environment* env)const;
		
		void updateEdgeFrequency();
		void selectSolutionsFromBasin(const HNSWbasin& model, Environment* env, Random* rnd);



		void generateSolutionsInBasinCross(
			const std::vector<int>& shuffleIds,
			std::vector<std::shared_ptr<SolutionType>>& sonSols,
			std::vector<bool>& activeFlag,
			ThreadInfo& thrdInfo,
			std::mutex& mtx,
			int& maxSolId,
			const HNSWbasin& model, Environment* env);


		void generateRandomSolutions(const HNSWbasin& model,
			Environment* env, Random* rnd);
		void generateRandomSolutionOneThread(
			ThreadInfo& thrdInfo,
			std::mutex& mtx,
			int & maxLoop,
			const HNSWbasin& model, Environment* env);
		void generateRandomSolutionMultiThread(
			const HNSWbasin& model, Environment* env, Random* rnd);

		void judgeInsideBasinOneThread(
			std::vector<SolutionType*>& sols, std::vector<bool>& flagInside,
			std::mutex& mtx, int& maxSolId, 
			const HNSWbasin& model, Environment* env, Random* rnd
		);
		void judgeInsideBasinMultiThread(
			std::vector<SolutionType*>& sols, std::vector<bool>& flagInside,
			const HNSWbasin& model, Environment* env, Random* rnd
		);
		
	public:
		void setBasinId(int basinId) {
			m_curBasinId = basinId;
		}
		const std::vector<SolutionType*>& curPop() { return m_curPop; }

		void initializeSingleThread(const HNSWbasin& model, Environment* env, Random* rnd) ;
		void initializeMultiThread(const HNSWbasin& model, Environment* env, Random* rnd);

		void evolveSingleThread(const HNSWbasin& model, Environment* env, Random* rnd);
		void evolveMultiThread(const HNSWbasin& model, Environment* env, Random* rnd);

		void setNumberThreads(int numThread, Random* rnd);

		
		bool terminationCondition() {
			return m_terminatedFlag;
		}

		SolutionType* bestSolution()const {
			return m_best;
		}
		void updateBestSolution() {
			if (m_curPop.empty())return;
			m_best = m_curPop.front();
			for (auto& it : m_curPop) {
				if (m_best->fitness() < it->fitness()) {
					m_best = it;
				}
			}
		}
		const std::vector<SolutionType*>& curPop()const {
			return m_curPop;
		}

	};
}

#endif // !OFEC_EAX_TSP_BASIN_POP_H

