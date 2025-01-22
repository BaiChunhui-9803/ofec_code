#include "nbn_improve_accuracy_multithread.h"
#include <shared_mutex>
#include <mutex>
#include <memory>

namespace ofec {
	namespace nbn {
	
		struct GlobalInfo {
			std::mutex m_globalMtx;
			std::vector<int> m_sortedIds;
			int m_curIdFrom = 0;
			int m_curIdTo = 0;
			bool m_finish = 0;
			int m_runningThread= 0;
			
			
		};



		void updateNBNinfoThreadTask(int curSolId, 
			NBNinfo& info,
			ofec::Environment* env, ofec::Random* rnd) {
			auto& nbnInfo = info;
		//	auto& updateIds = info.m_updateIds;
			auto curSol = nbnInfo.m_solbases[curSolId];

			
			for (int id = 0; id < curSolId; ++id) {
				auto otherSol = nbnInfo.m_solbases[id];
				/*if (curSol->fitness() > otherSol->fitness())*/ {
					double dis = curSol->variableDistance(otherSol->variableBase(), env);
					{
						//info.m_mtx.lock();

						if (nbnInfo.m_dis2parent[id] > dis) {
							//updateIds[id] = true;
							nbnInfo.m_dis2parent[id] = dis;
							nbnInfo.m_belong[id] = curSolId;

						}
	/*					else if (nbnInfo.m_dis2parent[id] == dis && rnd->uniform.next() < 0.5) {
							nbnInfo.m_dis2parent[id] = dis;
							nbnInfo.m_belong[id] = curSolId;
						}*/
						//info.m_mtx.unlock();
						
					}
				}
				//else if (curSol->fitness() < otherSol->fitness()) {

				//}
			}
		}


		void updateNBNinfoSingleThread(
			GlobalInfo& globalInfo,
			NBNinfo& info,
			ofec::Environment* env, ofec::Random* rnd
		) {

			//auto curRnd = std::make_shared<ofec::Random>(rnd->uniform.next());

			int curTaskId = 0;
			while (true) {
				
				{
					std::lock_guard<std::mutex> guard(globalInfo.m_globalMtx);
					if (globalInfo.m_finish)
					{
						--globalInfo.m_runningThread;
						curTaskId = -1;
					}
					else {
						curTaskId = globalInfo.m_curIdFrom;
						globalInfo.m_curIdFrom++;
						//globalInfo.m_totalTask.pop_back();
						
					}
				}
				if (curTaskId >= 0&&curTaskId< globalInfo.m_curIdTo) {
					updateNBNinfoThreadTask(curTaskId, info, env, rnd);
				}
				else break;
			}
		}


		void mergerInfoThreadTask(
			NBNinfo& curNbnInfo,
			std::vector<bool>& updateIds,
			std::vector<NBNinfo>& threadNBNinfos,
			int from, int to,
			ofec::Random*rnd
			) {
			for (int idx(from); idx < to; ++idx) {

				for (auto& oNbnInfo : threadNBNinfos) {
					if (curNbnInfo.m_dis2parent[idx] > oNbnInfo.m_dis2parent[idx]) {
						curNbnInfo.m_dis2parent[idx] = oNbnInfo.m_dis2parent[idx];
						curNbnInfo.m_belong[idx] = oNbnInfo.m_belong[idx];
						updateIds[idx] = true;
					}
					else if (curNbnInfo.m_dis2parent[idx] == oNbnInfo.m_dis2parent[idx] && rnd->uniform.next() < 0.5) {
						curNbnInfo.m_dis2parent[idx] = oNbnInfo.m_dis2parent[idx];
						curNbnInfo.m_belong[idx] = oNbnInfo.m_belong[idx];
					}
				}

				//= updateIds[idx]|| oUpdateIds[idx];
			}
		}

		void mergerInfoMultiThread(
			NBNinfo& curNbnInfo,
			std::vector<bool>& updateIds,
			std::vector<NBNinfo>& threadNBNinfos,
			std::vector<std::shared_ptr<ofec::Random>>& rnds) {

			int num_task = threadNBNinfos.size();
			int num_samples = curNbnInfo.m_belong.size();
			std::vector<std::thread> thrds;
			std::vector<int> tasks;
			UTILITY::assignThreads(num_samples, num_task, tasks);
			std::pair<int, int> from_to;


			for (size_t i = 0; i < num_task; ++i) {
				from_to.first = tasks[i];
				from_to.second = tasks[i + 1];
				thrds.push_back(std::thread(
					mergerInfoThreadTask,
					std::ref(curNbnInfo), std::ref(updateIds), std::ref(threadNBNinfos), from_to.first, from_to.second, rnds[i].get()));
			}

			for (auto& thrd : thrds)
				thrd.join();
			
		}


		


		void udpateNBNinfoIterationMultiThread(
			NBNinfo& m_nbnInfo,
			std::vector<bool>& updateIds,
			int from,int to,
			const std::vector<int>& sortedIds,
			ofec::Environment* env, ofec::Random* rnd) {

			updateIds.resize(m_nbnInfo.m_belong.size());
			std::fill(updateIds.begin(), updateIds.end(), false);
			
			int numThread = std::thread::hardware_concurrency();
			std::vector<NBNinfo> threadNBNinfos(numThread);
			std::vector<std::shared_ptr<ofec::Random>> rnds(numThread);
			for (auto& it : threadNBNinfos) {
				it.copy(m_nbnInfo);
			}
			for (auto& it : rnds) {
				it.reset(new Random(rnd->uniform.next()));
			}


			{
				GlobalInfo globalInfo;
				globalInfo.m_curIdFrom = from;
				globalInfo.m_curIdTo = to;
				globalInfo.m_sortedIds = sortedIds;
				globalInfo.m_finish = false;
				globalInfo.m_runningThread = numThread;
				std::vector<std::thread> thrds;

				for (size_t i = 0; i < numThread; ++i) {
					thrds.push_back(std::thread(
						updateNBNinfoSingleThread,
						std::ref(globalInfo), std::ref(threadNBNinfos[i]), env, rnds[i].get()));
				}

				for (auto& thrd : thrds)
					thrd.join();
				
			}

			mergerInfoMultiThread(m_nbnInfo,updateIds,threadNBNinfos,rnds);
			
			
		}

		
		//void udpate

		//void updateNBNinfo(int cursolId, 
	}
}