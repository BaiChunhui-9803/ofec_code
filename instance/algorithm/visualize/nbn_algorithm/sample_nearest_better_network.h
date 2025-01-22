#ifndef SAMPLE_NEAREST_BETTER_NETWORK_H
#define SAMPLE_NEAREST_BETTER_NETWORK_H

#include "../../../core/problem/solution.h"
#include "../../../utility/random/newran.h"
#include<memory>


namespace ofec {

	using NearestBetterNetworkSol = Solution<>;

	struct NearestBetterNetworkNode;
	struct NearestBetterNetworkListNode {
		NearestBetterNetworkListNode* m_before = nullptr;
		NearestBetterNetworkListNode* m_after = nullptr;
		std::shared_ptr<NearestBetterNetworkNode> m_info = nullptr;

		void clear() {
			m_before = nullptr;
			m_after = nullptr;
			m_info = nullptr;
		}

		void setInfo(const std::shared_ptr<NearestBetterNetworkNode>& cur) {
			m_info = cur;
		}
	};


	struct NearestBetterNetworkList {
		NearestBetterNetworkListNode m_head;
		NearestBetterNetworkListNode m_end;
		NearestBetterNetworkList() {
			clear();
		}
		void clear() {
			m_head.m_after = &m_end;
			m_end.m_before = &m_head;
		}
		void insertNode(NearestBetterNetworkListNode* cur) {
			insertNode(&m_head, cur);
			//cur->m_after = m_head.m_after;
			//m_head.m_after->m_before = cur;
			//cur->m_before = &m_head;
			//m_head.m_after = cur;
		}


		bool isEnd(NearestBetterNetworkListNode* cur) {
			return cur == &m_end;
		}

		NearestBetterNetworkListNode* head() {
			return &m_head;
		}
		NearestBetterNetworkListNode* front() {
			return m_head.m_after;
		}
		bool isHead(NearestBetterNetworkListNode* cur) {
			return cur == &m_head;
		}
		NearestBetterNetworkListNode* back() {
			return m_end.m_before;
		}
		
		


		static void deleteNode(NearestBetterNetworkListNode* cur) {
			cur->m_before->m_after = cur->m_after;
			cur->m_after->m_before = cur->m_before;
			cur->m_after = cur->m_before = nullptr;
		}


		static void insertNode(NearestBetterNetworkListNode* ahead,NearestBetterNetworkListNode* cur) {
			cur->m_after = ahead->m_after;
			ahead->m_after->m_before = cur;
			cur->m_before = ahead;
			ahead->m_after = cur;
		}


	};



	class NBNStaticPriorityQueue {
	public:
		using Tdata = std::shared_ptr<NearestBetterNetworkNode>;
	protected:

		int m_problem.get() = -1;
		int m_random.get() = -1;
		Tdata m_cur_sol;
		int m_cur_size = 0;
		const int m_max_size = 50;
		std::array<Tdata, 51> m_data;
		std::function<bool(
			const Tdata& cur,
			const  Tdata& a, const Tdata& b,
			Problem *pro, Random *rnd)> m_com_fun;
	public:



		double avgDis();

		int size() const{
			return m_cur_size;
		}

		const std::array<Tdata, 51>& data()const {
			return m_data;
		}
		void clearData(const Tdata& data) {
			std::fill(m_data.begin(), m_data.end(), data);
			m_cur_size = 0;
			m_cur_sol = data;
		}

		void setFunction(const std::function<bool(
			const Tdata& cur,
			const  Tdata& a, const Tdata& b,Problem *pro, Random *rnd)> &com_fun) {
			m_com_fun = com_fun;
		}
		void setCurSol(const Tdata& cur) {
			m_cur_sol = cur;
		}
		void setIdPro(Problem *pro) {
			m_problem.get() = pro;
		}
		void setIdRnd(Random *rnd) {
			m_random.get() = rnd;
		}
		bool empty() {
			return m_cur_size == 0;
		}
		void clear() {
			m_cur_size = 0;
		}
		bool push(const Tdata& a, Tdata& b) {
			m_data[m_cur_size++] = a;
			int cur_size = (m_cur_size - 1);
			while (cur_size) {
				if (m_com_fun(
					m_cur_sol,
					m_data[cur_size], m_data[(cur_size - 1) / 2.0],
					m_problem.get(), m_random.get()
					)) {
					swap(m_data[cur_size], m_data[(cur_size - 1) / 2.0]);
					cur_size = (cur_size - 1) / 2.0;
				}
				else break;
			}
			if (m_cur_size == m_max_size+1) {
				b = m_data[m_cur_size-1];
				m_cur_size = m_max_size;
				return true;
			}
			else return false;
		}
		void pop() {
			--m_cur_size;
			swap(m_data[0], m_data[m_cur_size]);
			int cur_size(0);
			while (cur_size < m_cur_size) {
				int left = cur_size * 2 + 1;
				int right = cur_size * 2 + 2;
				if (left >= m_cur_size) break;
				if (right >= m_cur_size) {
					if (m_com_fun(
						m_cur_sol,
						m_data[left], m_data[cur_size],
						m_problem.get(), m_random.get()
						)) {
						swap(m_data[left], m_data[cur_size]);
						m_cur_size = left;
					}
					break;
				}
				else {
					int large_idx(left);
					if (m_com_fun(
						m_cur_sol,
						m_data[right], m_data[left],
						m_problem.get(), m_random.get()
						)) {
						large_idx = right;
					}
					if (m_com_fun(
						m_cur_sol,
						m_data[large_idx], m_data[cur_size],
						m_problem.get(), m_random.get()
						)) {
						swap(m_data[large_idx], m_data[cur_size]);
						cur_size = large_idx;
					}
					else break;
				}
			}
		}
		const Tdata& top()const {
			return m_data[0];
		}
	};




	struct  NearestBetterNetworkNode {
		using SolutionType = Solution<>;
		using Tdata = std::shared_ptr<NearestBetterNetworkNode>;
		int m_sample_id = -1;
		int m_sorted_id = -1;
		double m_disToParent = std::numeric_limits<double>::max();
		std::shared_ptr<NearestBetterNetworkNode> m_parent = nullptr;
		NearestBetterNetworkListNode m_parent_list_node;
		NearestBetterNetworkList m_sons;
		NearestBetterNetworkListNode m_sorted_list_node;
		std::shared_ptr<SolutionType> m_cur_sol = nullptr;

		NBNStaticPriorityQueue m_neighbors;
		double m_neighbors_sum_dis = 0;
		unsigned long long m_visited_id = 0;
		bool m_flag_opt = false;

		void clear() {
			m_sample_id = -1;
			m_sorted_id = -1;
			m_disToParent = std::numeric_limits<double>::max();

			m_parent = nullptr;
			m_parent_list_node.clear();
			m_sons.clear();
			m_sorted_list_node.clear();
			m_cur_sol = nullptr;
			m_neighbors.clearData(nullptr);
			m_neighbors_sum_dis = 0;
			m_visited_id = 0;
		}

		bool flagNoParent() {
			return m_parent == nullptr;
		}

		bool updateParent(std::shared_ptr<NearestBetterNetworkNode>& parent,Problem *pro) {

			if (parent->m_cur_sol->fitness() > m_cur_sol->fitness()) {
				double dis(m_cur_sol->variableDistance(*parent->m_cur_sol, pro));
				if (dis < m_disToParent) {
					if (m_parent != nullptr) {
						NearestBetterNetworkList::deleteNode(&m_parent_list_node);
					}
					m_disToParent = dis;
					m_parent = parent;
					m_parent->m_sons.insertNode(&m_parent_list_node);
					return true;
				}
			}
			return false;
		}

		void setSol(std::unique_ptr<SolutionType>& sol) {
			m_cur_sol = std::move(sol);
		}
		void setSampleId(int id) {
			m_sample_id = id;
		}
		void initialize(
			Problem *pro,
			Random *rnd,
			int id,
			std::shared_ptr<SolutionType>& sol,
			const std::function<bool(
				const Tdata &cur, 
				const Tdata& a, const Tdata& b, 
				Problem *pro, Random *rnd)>& com_fun,
			const std::shared_ptr<NearestBetterNetworkNode>& curNode
			) {

			m_sample_id = id;
			m_sorted_id = -1;
			m_disToParent = std::numeric_limits<double>::max();

			m_parent = nullptr;
			m_parent_list_node.setInfo(curNode);
			m_sons.clear();
			m_sorted_list_node.setInfo(curNode);
			m_cur_sol = sol;
			m_neighbors.clearData(nullptr);
			m_neighbors.setIdPro(pro);
			m_neighbors.setIdRnd(rnd);
			m_neighbors.setFunction(com_fun);
			m_neighbors.setCurSol(curNode);
			
			m_neighbors_sum_dis = 0;
			m_visited_id = 0;
		}

		bool addNeighbor(Problem *pro,const std::shared_ptr<NearestBetterNetworkNode>& sol) {
			std::shared_ptr<NearestBetterNetworkNode> del_node = nullptr;
			m_neighbors_sum_dis += sol->m_cur_sol->variableDistance(*m_cur_sol, pro);
			if (m_neighbors.push(sol,del_node)) {
				m_neighbors_sum_dis -= del_node->m_cur_sol->variableDistance(*m_cur_sol,pro);
				return true;
			}
			return false;
		}
	};



	struct SampleNearestBetterNetworkRecord {
		using SolutionType = Solution<>;
		using Tdata = std::shared_ptr<NearestBetterNetworkNode>;
		std::vector<std::shared_ptr<SolutionType>> m_samples;
		std::vector<int> m_sorted_ids;
		std::vector<int> m_parent_ids;
		std::vector<std::vector<int>> m_neighbor_ids;
		//std::vector<double> m_fitness;
		int m_random.get();
		int m_problem.get();
		//std::function<void(SolutionType& sol,Problem *pro)> m_eval_fun;
		std::function<bool(
			const Tdata& cur,
			const  Tdata& a, const Tdata& b,
			Problem *pro, Random *rnd
			)> m_com_fun;
		//unsigned long long m_cur_visited_id = 0;
		void getNearestBetterNetwork(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong)
		{
			sols.resize(m_samples.size());
			fitness.resize(m_samples.size());
			belong.resize(m_samples.size());
			for (int idx(0); idx < sols.size(); ++idx) {
				sols[idx].reset(new SolutionType(*m_samples[idx]));
				fitness[idx] = sols[idx]->fitness();
				belong[idx] = m_parent_ids[idx];
			}
		}

		bool isEqualNetwork(const SampleNearestBetterNetworkRecord& other) const{

			if (m_samples.size() != other.m_samples.size()) {
				std::cout << "error: different_size" << std::endl;
				return false;
			}
			if (m_sorted_ids != other.m_sorted_ids) {
				std::cout << " sort idxs different size" << std::endl;
			}
			if (m_parent_ids != other.m_parent_ids) {
				std::cout << "parent ids different" << std::endl;
			}

			std::vector<int> vi1;
			std::vector<int> vi2;
			for (int idx(0); idx < m_samples.size(); ++idx) {
				vi1 = m_neighbor_ids[idx];
				vi2 = other.m_neighbor_ids[idx];
				std::sort(vi1.begin(), vi1.end());
				std::sort(vi2.begin(), vi2.end());
				if (vi1 != vi2) {
					std::cout << "different neighbor idx \t" << idx << "\t" << std::endl;
				}
			}


		}
	};

	
	class SampleNearestBetterNetwork {
	public :
		using SolutionType = Solution<>;
		using Tdata = std::shared_ptr<NearestBetterNetworkNode>;
	protected:

		int m_random.get();
		int m_problem.get();
		std::function<bool(
			const Tdata& cur,
			const  Tdata& a, const Tdata& b,
			Problem *pro, Random *rnd
			)> m_com_fun;
	
		
		std::vector<std::shared_ptr<NearestBetterNetworkNode>> m_samples;
		NearestBetterNetworkList m_sorted_list;
		unsigned long long m_cur_visited_id = 0;

	protected:
		std::shared_ptr<NearestBetterNetworkNode>& insertNode(std::shared_ptr<SolutionType>& sol) {
			int sampleId(m_samples.size());
			m_samples.emplace_back(new NearestBetterNetworkNode);
			m_samples.back()->initialize(m_problem.get(), m_random.get(),
				sampleId, sol, 
				m_com_fun,m_samples.back());
			m_samples.back()->m_cur_sol->setId(sampleId);
			return m_samples.back();
		}
		void setComFun() {
			m_com_fun = []
			(const Tdata& cur, const Tdata& a, const Tdata& b,
				Problem *pro, Random *rnd) {
				double disA = cur->m_cur_sol->
					variableDistance(*a->m_cur_sol, pro);
				double disB = cur->m_cur_sol->
					variableDistance(*a->m_cur_sol, pro);
				if (disA == disB) {
					if (rnd->uniform.next() < 0.5) return true;
					else return false;
				}
				else return disA < disB;
			};
		}
	public:

		int idRnd()const { return m_random.get();}
		int idPro()const { return m_problem.get();}
		const std::function<bool(
			const Tdata& cur,
			const  Tdata& a, const Tdata& b,
			Problem *pro, Random *rnd
			)>& comFun()const {
			return m_com_fun;
		}

		void clear() {
			for (auto& it : m_samples) {
				it->clear();
			}
			m_sorted_list.clear();
			m_cur_visited_id = 0;
			m_samples.clear();
		}
	
		const std::vector<std::shared_ptr<NearestBetterNetworkNode>>& getSamples() {
			return m_samples;
		}

		//void setEvalFun(const std::function<void(SolutionType& sol, Problem *pro)>& eval_fun) {
		//	m_eval_fun = eval_fun;
		//}
		~SampleNearestBetterNetwork() {
			clear();
		}
		void initialize(Random *rnd, Problem *pro) {
			clear();
			m_random.get() = rnd;
			m_problem.get() = pro;
			m_cur_visited_id = 0;
			setComFun();
			
		}
		bool addRandomSol(std::shared_ptr<SolutionType>& sol,bool flag_opt=false);
		void addNearSol(const std::shared_ptr<NearestBetterNetworkNode>& center, std::unique_ptr<SolutionType>& sol);
		void initRandomSol(int initSols, const std::function<void(SolutionType& sol, Problem *pro)>& eval_fun);
		void sampleRandom(
			std::vector<std::shared_ptr<SolutionType>>& sols, const std::function<void(SolutionType& sol, Problem *pro)>& eval_fun);
		void sampleRandomThreadTask(
			std::vector<std::shared_ptr<SolutionType>>& sols,
			std::pair<int, int> from_idx,double randSeed);
		void updateSortedId() {
			auto curIter(m_sorted_list.front());
			int cur_id(0);
			while (!m_sorted_list.isEnd(curIter)) {
				curIter->m_info->m_sorted_id = cur_id++;
				curIter = curIter->m_after;
			}
		}

		void getNearestBetterNetwork(
			std::vector<std::shared_ptr<SolutionType>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong);



		void getNearestBetterNetwork(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong);

		void getNearestBetterNetworkShareMemory(
			std::vector<std::shared_ptr<SolutionType>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong);

		void outputData(SampleNearestBetterNetworkRecord& record);
		void insertData(const SampleNearestBetterNetworkRecord& record);
	};

}



#endif
