#ifndef NBN_NEAREST_BETTER_CALCULATOR_MULTITHREAD_V2_H
#define NBN_NEAREST_BETTER_CALCULATOR_MULTITHREAD_V2_H

#include <vector>
#include <memory>
#include <queue>
#include <mutex>
#include <thread>
#include <algorithm>
#include <shared_mutex>
#include <chrono>


#include "../../../core/problem/solution.h"
#include "../../../core/problem/problem.h"
#include "../../../utility/function/custom_function.h"
#include "../../../utility/function/visit_lazy_list.h"

namespace ofec {
	struct NBN_NearestBetterCalculatorMultiThread_v2 {

	protected: 


		static bool compareLarger(int curIdA, double value1, int curIdB, double value2) {
			if (value1 == value2) {
				return curIdA > curIdB;
			}
			else return value1 > value2;
		}

		struct NodeValueStruct {
			int m_cur_id = 0;
			double value = 0;
			bool operator<(const NodeValueStruct& a) const
			{
				return compareLarger(m_cur_id, value, a.m_cur_id, a.value);
				//return value > a.value;
			}
			NodeValueStruct() = default;
			NodeValueStruct(int id, double curvalue):m_cur_id(id), value(curvalue){}
			void setId(int id) {
				m_cur_id = id;
			}

		};

		struct NodeInfo {

		private:
			std::shared_mutex m_after_mutex;
			bool m_finished = false;
			std::vector<NodeValueStruct> m_after_nodes;
		public:
			
			int m_stage = 0;
			std::priority_queue<NodeValueStruct>* nei_que = nullptr;
			std::vector<int>* attraction = nullptr;
			custom_fun::VisitedLazyList* visited = nullptr;

			std::vector<int> m_betterRange;

			bool isFinish() {
				std::shared_lock lock(m_after_mutex);
				return m_finished;
			}
			
			void insertAfterNode(const NodeValueStruct& cur) {
				std::unique_lock lock(m_after_mutex);
				m_after_nodes.push_back(cur);
			}

			void finishTask(
				std::vector<NodeValueStruct>& after_nodes, 
				std::vector<int>& betterRange) {
				std::unique_lock lock(m_after_mutex);
				after_nodes = std::move(m_after_nodes);
				m_betterRange = std::move(betterRange);
				m_finished = true;
			}
			

		};
		
		struct GlobalInfo {

		private:
			//std::shared_mutex m_belongMutex;
			std::vector<int>& m_belong;
			std::vector<double>& m_dis2parent;
			
			std::shared_mutex m_mutex;
			int m_leftTask = 0;
			std::priority_queue<NodeValueStruct> m_tasks;

			std::mutex m_visitedQueMtx;
			std::queue<custom_fun::VisitedLazyList*> visitedQue;
		public:


			std::mutex m_memory_mtx;

			const std::vector<SolutionBase*>& m_sols;
			const std::vector<double>& m_fitness;
			const std::vector<std::vector<int>>& m_neighbors;

			std::vector<NodeInfo*> nodeInfos;



			GlobalInfo(
				const std::vector<SolutionBase*>& sols,
				const std::vector<double>& fitness,
				const std::vector<std::vector<int>>& neighbors,
				std::vector<int>& belong,
				std::vector<double>& dis2parent) :
				m_sols(sols), m_fitness(fitness),
				m_neighbors(neighbors), m_belong(belong), m_dis2parent(dis2parent){}

			bool compareTwoSolsLarger(int a, int b) {
				return compareLarger(a, m_fitness[a], b, m_fitness[b]);
			}
			int getBelong(int solId) {
				//std::shared_lock lock(m_belongMutex);
				return m_belong[solId];
			}
			void setBelongDistance(int solId, int belong, double distance) {
				//std::unique_lock lock(m_belongMutex);
				m_belong[solId] = belong;
				m_dis2parent[solId] = distance;
			}

			void addTaskSingleTask(const NodeValueStruct& curTask) {
				m_tasks.push(curTask);
				++m_leftTask;
			}

			void removeTask() {
				std::shared_lock lock(m_mutex);
				--m_leftTask;
			}
			int leftTask() {
				std::shared_lock lock(m_mutex);
				return m_leftTask;
			}

			void getTask(NodeValueStruct& curTask) {
				std::unique_lock lock(m_mutex);
				if (!m_tasks.empty()) {
					curTask = m_tasks.top();
					m_tasks.pop();
				}
			}
			void insertTask(std::vector<NodeValueStruct>& tasks) {
				std::unique_lock lock(m_mutex);
				while (!tasks.empty()) {
					m_tasks.emplace(std::move(tasks.back()));
					tasks.pop_back();
				}
			}
			void insertTask(NodeValueStruct& task) {
				std::unique_lock lock(m_mutex);
				m_tasks.emplace(std::move(task));
			}

			void finishTask() {
				std::unique_lock lock(m_mutex);
				--m_leftTask;
			}

			void getVisitedQue(custom_fun::VisitedLazyList*& visited) {
				m_visitedQueMtx.lock();
				if (!visitedQue.empty()) {
					visited = visitedQue.front();
					visitedQue.pop();
				}
				m_visitedQueMtx.unlock();
			}

			void insertVisitedQue(custom_fun::VisitedLazyList*& visited) {
				m_visitedQueMtx.lock();
				if (!visitedQue.empty()) {
					visitedQue.push(visited);
					visited = nullptr;
					visitedQue.pop();
				}
				m_visitedQueMtx.unlock();
			}

			void assignMemory(int size) {
				nodeInfos.resize(size);
				for (auto& it : nodeInfos) {
					it = new NodeInfo;
				}
			}


			void releaseMemory() {
				for (auto& it : nodeInfos) {
					delete it;
				}
				nodeInfos.clear();
				while (!visitedQue.empty()) {
					delete visitedQue.front();
					visitedQue.pop();
				}
			}

		};


	protected:



		static void resetVisited(std::vector<unsigned long long>& visited,
			unsigned long long& curVisited) {
			curVisited = 0;
			std::fill(visited.begin(), visited.end(), curVisited);
		}


		static void updateParentTask(GlobalInfo& globalInfo,
			ofec::Problem* pro) {
			NodeValueStruct curTask;
			custom_fun::VisitedLazyList* visitedList = nullptr;
			bool finishTask = false;
			while (true) {

				
				if (globalInfo.leftTask() == 0)break;
				curTask.setId(-1);
				globalInfo.getTask(curTask);
				if (curTask.m_cur_id == -1) {
					std::this_thread::sleep_for(std::chrono::seconds(2));
				}
				else {
					finishTask = false;
					auto& curInfo = *globalInfo.nodeInfos[curTask.m_cur_id];
					auto& curSolId = curTask.m_cur_id;

					// allocate memory
					if (curInfo.m_stage == 0) {

						bool flag = false;
						globalInfo.m_memory_mtx.lock();
						{
							int number = 2;
							if (curInfo.visited == nullptr) {
								++number;
							}
							number *= globalInfo.nodeInfos.size();
							auto buffer = std::get_temporary_buffer<NodeValueStruct>(number);
							if (buffer.second >= number) {
								flag = true;
							}
						}
						if (flag) {
							curInfo.nei_que = new std::priority_queue<NodeValueStruct>();
							curInfo.attraction = new std::vector<int>();
							globalInfo.getVisitedQue(curInfo.visited);

							if (curInfo.visited == nullptr) {
								curInfo.visited = new custom_fun::VisitedLazyList;
								curInfo.visited->resize(globalInfo.nodeInfos.size());
							}
						}
						globalInfo.m_memory_mtx.unlock();
						if (!flag) {
							globalInfo.insertTask(curTask);
							//break;
						}
						else {
							curInfo.m_stage = 1;
						}
					}
			
					if (curInfo.m_stage == 1) {

						auto& nei_que = *curInfo.nei_que;
						auto& cur_sol(globalInfo.m_sols[curSolId]);
						NodeValueStruct curQueNode;
						curInfo.visited->Reset();
						curInfo.visited->MarkAsVisited(curSolId);

						auto curNeibor = globalInfo.m_neighbors[curSolId];

						for (auto& nei_info : curNeibor) {
							auto& nei_idx = nei_info;
							if (curInfo.visited->NotVisited(nei_idx)) {
								curInfo.visited->MarkAsVisited(nei_idx);
								curQueNode.m_cur_id = nei_idx;
								{
									auto& nei_sol(globalInfo.m_sols[nei_idx]);
									curQueNode.value = pro->normalizedVariableDistance(*cur_sol, *nei_sol);
									//cur_node.m_dis2parent = cur_sol->norVariableDistance(*par_sol, m_id_pro);
								}
								nei_que.push(curQueNode);
							}
						}
						curInfo.m_stage = 2;
					}

					if (curInfo.m_stage == 2) {
						auto& nei_que = *curInfo.nei_que;
						auto& cur_sol(globalInfo.m_sols[curSolId]);
						NodeValueStruct curQueNode;
						auto& attraction = *curInfo.attraction;
						int nei_node = -1;
						while (!nei_que.empty()) {
							curQueNode = nei_que.top();
							nei_node = curQueNode.m_cur_id;
							
							if (globalInfo.compareTwoSolsLarger(nei_node, curSolId)) {
								curInfo.m_stage = 3;
								break;
							}

							std::vector<int> betterRange;
							

							
							if (!globalInfo.nodeInfos[nei_node]->isFinish()) {
								betterRange = globalInfo.m_neighbors[nei_node];
							}
							else {
								betterRange = globalInfo.nodeInfos[nei_node]->m_betterRange;
							}

							if (globalInfo.getBelong(nei_node) != -1){
								nei_que.pop();
								attraction.push_back(curQueNode.m_cur_id);
								for (auto& son_idx : betterRange) {
									if (curInfo.visited->NotVisited(son_idx)) {
										curInfo.visited->MarkAsVisited(son_idx);
										curQueNode.m_cur_id = son_idx;
										auto& son_sol(globalInfo.m_sols[son_idx]);
										curQueNode.value = pro->normalizedVariableDistance(*cur_sol, *son_sol);
										nei_que.push(curQueNode);
									}
								}
							}
							else {
								curInfo.m_stage = 3;
								break;
							}
						}
						if (curInfo.m_stage == 2) {
							globalInfo.nodeInfos[nei_node]->insertAfterNode(curTask);
						}
					}

					if (curInfo.m_stage == 3) {

						auto& nei_que = *curInfo.nei_que;
						auto& cur_sol(globalInfo.m_sols[curSolId]);
						NodeValueStruct curQueNode;
						std::vector<int> betterRange;
						if (!nei_que.empty()) {
							NodeValueStruct parentNode = nei_que.top();
							while (!nei_que.empty()) {
								curQueNode = nei_que.top();
								nei_que.pop();
								betterRange.push_back(curQueNode.m_cur_id);
							}
							if (parentNode.m_cur_id != -1) {
								int parSolId = parentNode.m_cur_id;
								globalInfo.setBelongDistance(curSolId, parSolId, parentNode.value);
							}
						}
						std::vector<NodeValueStruct> afterNodes;
						{
							globalInfo.nodeInfos[curTask.m_cur_id]->finishTask(afterNodes, betterRange);
						}
						globalInfo.insertTask(afterNodes);


						{
							delete curInfo.nei_que;
							curInfo.nei_que = nullptr;
							delete curInfo.attraction;
							curInfo.attraction = nullptr;

							globalInfo.insertVisitedQue(curInfo.visited);
						}

						globalInfo.removeTask();
					}
				}
			}
		}


	public:

		static void calculate(const std::vector<SolutionBase*>& m_sols,
			const std::vector<double>& fitness,
			const std::vector<std::vector<int>>& neighbors,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			ofec::Problem* pro, ofec::Random* rnd) {
			std::vector<unsigned long long> visited(m_sols.size(), 0);
			unsigned long long cur_visited = 0;
			std::vector<std::vector<int>> better_range_nodeId(m_sols.size());
			belong.resize(m_sols.size());
			std::fill(belong.begin(), belong.end(), -1);

			dis2parent.resize(m_sols.size());
			std::fill(dis2parent.begin(), dis2parent.end(), std::numeric_limits<double>::max());

			std::vector<int> sorted_idx(m_sols.size());
			for (int idx(0); idx < sorted_idx.size(); ++idx) {
				sorted_idx[idx] = idx;
			}
			//	auto& nbn_nodes(curTreeNode.m_division_nodes);
			std::sort(sorted_idx.begin(), sorted_idx.end(), [&](
				int a, int  b
				) {
				return fitness[a] < fitness[b];
			});

			//curTreeNode.m_cur_visited2 = 0;
			//for (auto& curNode : curTreeNode.m_division_nodes) {
			//	curNode.m_better_range_nodeId.clear();
			//	curNode.m_better_range_nodeId.push_back(curNode.m_node_id);
			//}
			for (int sidx(0); sidx + 1 < sorted_idx.size(); ++sidx) {
				//std::cout << "curId\t" << sidx << std::endl;
				auto sortedId = sorted_idx[sidx];
				//	auto& cur_node(nbn_nodes[sortedId]);

				updateParent(m_sols, fitness, neighbors, sortedId,
					better_range_nodeId, belong, dis2parent, visited, cur_visited,
					pro, rnd);

			}

			belong[sorted_idx.back()] = sorted_idx.back();

			//	curTreeNode.m_best_divNodeId = sorted_idx.back();
			//	curTreeNode.m_best_solId = nbn_nodes[sorted_idx.back()].m_representative_solId;

			for (int idx(0); idx < belong.size(); ++idx) {
				if (belong[idx] == -1) belong[idx] = idx;
			}
		}

		



		static void calculateMutiThread(const std::vector<SolutionBase*>& sols,
			const std::vector<double>& fitness,
			const std::vector<std::vector<int>>& neighbors,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			ofec::Problem* pro) {

			belong.resize(sols.size());
			std::fill(belong.begin(), belong.end(), -1);

			dis2parent.resize(sols.size());
			std::fill(dis2parent.begin(), dis2parent.end(), std::numeric_limits<double>::max());



			GlobalInfo globalInfo(sols,fitness,neighbors,belong,dis2parent);
			globalInfo.assignMemory(sols.size());
			NodeValueStruct curTask;

			int bestId = 0;
			{
				//double bestFit = fitness.front();
				for (int idx(0); idx < sols.size(); ++idx) {
					if (globalInfo.compareTwoSolsLarger(idx, bestId)) {
						bestId = idx;
					}
					//if (bestFit < fitness[idx]) {
					//	bestId = idx;
					//	bestFit = fitness[idx];
					//}
				}
			}

			for (int idx(0); idx < sols.size(); ++idx) {
				if (idx != bestId) {
					curTask.m_cur_id = idx;
					curTask.value = fitness[idx];
					globalInfo.addTaskSingleTask(curTask);
				}
			}
			int num_task = std::thread::hardware_concurrency();

			//num_task = 1;
			std::vector<std::thread> thrds;
			for (size_t i = 0; i < num_task; ++i) {
				thrds.push_back(std::thread(
					&NBN_NearestBetterCalculator::updateParentTask, std::ref(globalInfo), pro));
			}
			for (auto& thrd : thrds)
				thrd.join();


			belong[bestId] = bestId;
			for (int idx(0); idx < belong.size(); ++idx) {
				if (belong[idx] == -1) belong[idx] = idx;
			}


			

		}


		static void updateParent(const std::vector<SolutionBase*>& m_sols,
		const std::vector<double>& fitness,
			const std::vector<std::vector<int>>& neighbors,
			int curSolId,
			std::vector<std::vector<int>>& better_range_nodeId,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<unsigned long long>& visited,
			unsigned long long& curVisited,
			ofec::Problem* pro, ofec::Random* rnd) {
			if (++curVisited == 0) {
				resetVisited(visited, curVisited);
			}
			visited[curSolId] = curVisited;

			auto& cur_sol(m_sols[curSolId]);

			std::priority_queue<NodeValueStruct> nei_que;
			NodeValueStruct curQueNode;
			auto curNeibor = neighbors[curSolId];
			
			rnd->uniform.shuffle(curNeibor.begin(), curNeibor.end());
			for (auto& nei_info : curNeibor) {
				auto& nei_idx = nei_info;
				//auto nei_node(&nbn_nodes[nei_idx]);
				if (visited[nei_idx] != curVisited) {
					visited[nei_idx] = curVisited;
					//curQueNode.m_cur = nei_node;
					curQueNode.m_cur_id = nei_idx;
					{
						auto& nei_sol(m_sols[nei_idx]);
						curQueNode.value = pro->normalizedVariableDistance(*cur_sol, *nei_sol);
						//cur_node.m_dis2parent = cur_sol->norVariableDistance(*par_sol, m_id_pro);
					}
					nei_que.push(curQueNode);
				}
			}


			std::vector<int> attraction;
			while (!nei_que.empty()) {
				curQueNode = nei_que.top();
				auto nei_node = curQueNode.m_cur_id;

				if (belong[nei_node] != -1)
					//if (m_belong[nei_node->m_belong] != nei_node->m_representative_solId)
					/*if (m_fitness[cur_node.m_representative_solId] > m_fitness[nei_node->m_representative_solId])*/ {
					nei_que.pop();
					attraction.push_back(curQueNode.m_cur_id);
					for (auto& son_idx : better_range_nodeId[nei_node]) {
					//	auto son_node(&nbn_nodes[son_idx]);
						if (visited[son_idx] != curVisited) {
							visited[son_idx] = curVisited;
							//curQueNode.m_cur = son_node;
							curQueNode.m_cur_id = son_idx;
							{

								auto& son_sol(m_sols[son_idx]);
								curQueNode.value = pro->normalizedVariableDistance(*cur_sol, *son_sol);
								//cur_node.m_dis2parent = cur_sol->norVariableDistance(*par_sol, m_id_pro);
							}
							nei_que.push(curQueNode);
						}
					}
				}
				else break;
			}
			if(!nei_que.empty()){
				//cur_node.m_better_range_nodeId = peak_ranges;
				better_range_nodeId[curSolId].clear();
				NodeValueStruct parentNode = nei_que.top();
				while (!nei_que.empty()) {
					curQueNode = nei_que.top();
					nei_que.pop();
					better_range_nodeId[curSolId].push_back(curQueNode.m_cur_id);
				}

				
			//	cur_node.m_direct_parent_nodeId = cur_node.m_parent_nodeId = parentNode.m_cur_id;
				if (parentNode.m_cur_id != -1) {
				//	int curSolId = cur_node.m_representative_solId
					int parSolId = parentNode.m_cur_id;
					dis2parent[curSolId] = parentNode.value;
					belong[curSolId] = parSolId;

				}
			}
		}
		


		static void modifySols(const std::vector<SolutionBase*>& m_sols,
			const std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			ofec::Problem* pro, ofec::Random* rnd) {

			std::vector<int> solIds(m_sols.size());
			for (int idx(0); idx < solIds.size(); ++idx) {
				solIds[idx] = idx;
			}
			std::sort(solIds.begin(), solIds.end(), [&](int a,int b) {
				return compareLarger(a, fitness[a], b, fitness[b]);
				//return fitness[a] < fitness[b];
			});


			for (int idx(0); idx < solIds.size(); ++idx) {
				auto curSolId = solIds[idx];
				auto& cursol = m_sols[curSolId];
				if (belong[idx] == idx) {
					for (int idy(idx + 1); idy < solIds.size(); ++idy) {
						auto otherSolId = solIds[idy];
						auto& otherSol = m_sols[otherSolId];
						double dis = pro->normalizedVariableDistance(*cursol, *otherSol);
						if (dis < dis2parent[curSolId]) {
							dis2parent[curSolId] = dis;
							belong[curSolId] = otherSolId;
						}
					}
			    }
			}
		}
	};

}


#endif