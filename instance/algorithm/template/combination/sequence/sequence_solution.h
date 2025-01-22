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
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* Solution for combinatorial problem that can be transformed to problem with edges
*
*********************************************************************************/


#ifndef SEQUENCE_SOLUTION_H
#define SEQUENCE_SOLUTION_H
#include<vector>
namespace ofec {

	class SequenceSolution
	{
	public:
		SequenceSolution() = default;
		~SequenceSolution() = default;
		SequenceSolution& operator =(const SequenceSolution& rhs) = default;
		const std::vector<std::pair<int, int>>& currentEdges()const {
			return m_cur_edges;
		}
		std::vector<std::pair<int, int>>& currentEdges() {
			return m_cur_edges;
		}
		bool& flagCurrentEdgesUpdated() {
			return m_flag_cur_edges_updated;
		}
		bool flagCurrentEdgesUpdated() const{
			return m_flag_cur_edges_updated;
		}
		
		const std::vector<std::pair<int, int>>& edges()const {
			return m_edges;
		}
		std::vector<std::pair<int, int>>& edges() {
			return m_edges;
		}

		bool& flagEdgesUpdated() {
			return m_flag_edges_updated;
		}
		bool flagEdgesUpdated() const {
			return m_flag_edges_updated;
		}

		int getLoc()const {
			return m_loc;
		}
		int& getLoc() {
			return m_loc;
		}
		const std::vector<int>& feasibleNeighbors()const {
			return m_feasible_neighbors;
		}
		std::vector<int>& feasibleNeighbors() {
			return m_feasible_neighbors;
		}
		bool flagFeasibleNeighborsUpdated()const {
			return m_flag_feasible_neighbors_updated;
		}
		bool& flagFeasibleNeighborsUpdated() {
			return m_flag_feasible_neighbors_updated;
		}

		void setLocalSearch(bool flag) {
			m_flag_local_search = flag;
		}
		bool flagLocalSearch()const {
			return m_flag_local_search;
		}

		size_t& evaluateTimes() {
			return m_evaluate_times;
		}

		virtual void reset() {
			m_flag_local_search = false;
			m_feasible_neighbors.clear();
			m_flag_feasible_neighbors_updated = false;
			m_edges.clear();
			m_flag_edges_updated = false;
			m_cur_edges.clear();
			m_flag_cur_edges_updated = false;
			m_loc = -1;
			m_evaluate_times = 0;
		}

	//	virtual void updateInfo(Problem *pro, Algorithm *alg) {}
		
	protected:

		bool m_flag_local_search = false;
		std::vector<int> m_feasible_neighbors;
		bool m_flag_feasible_neighbors_updated = false;
		std::vector<std::pair<int, int>> m_edges;
		bool m_flag_edges_updated = false;
		std::vector<std::pair<int, int>> m_cur_edges;
		bool m_flag_cur_edges_updated = false;
		int m_loc = -1;
		size_t m_evaluate_times = 0;
	};


	template<typename Population>
	void calculateDensityMatrix(const Population& pop,
		std::vector<std::vector<ofec::Real>>& mat) {
		for (auto& it : mat) {
			for (auto& it2 : it) {
				it2 = 0;
			}
		}

		for (int idx(0); idx < pop.size(); ++idx) {
			for (auto& edge_iter : pop[idx].edges()) {
				++mat[edge_iter.first][edge_iter.second];
			}
		}
		ofec::Real max_val(0);
		for (auto& it : mat) {
			for (auto& it2 : it) {
				max_val = std::max(max_val, it2);
			}
		}

		if (max_val == 0) {
			int stop = -1;
		}
		for (auto& it : mat) {
			for (auto& it2 : it) {
				it2 /= max_val;
			}
		}
	}


	//template<typename TSolution>
	//void calculateDensityMatrix(const std::vector<std::unique_ptr<TSolution>>& indis,
	//	std::vector<std::vector<ofec::Real>>& mat) {
	//	for (auto& it : mat) {
	//		for (auto& it2 : it) {
	//			it2 = 0;
	//		}
	//	}
	//	for (auto& it : indis) {
	//		for (auto& edge_iter : it->edges()) {
	//			++mat[edge_iter.first][edge_iter.second];
	//		}
	//	}
	//	ofec::Real max_val(0);
	//	for (auto& it : mat) {
	//		for (auto& it2 : it) {
	//			max_val = std::max(max_val, it2);
	//		}
	//	}
	//	for (auto& it : mat) {
	//		for (auto& it2 : it) {
	//			it2 /= max_val;
	//		}
	//	}
	//}

}
#endif