
#ifndef DISTANCE_CALCULATOR_H
#define DISTANCE_CALCULATOR_H

#include"../../../../../core/definition.h"

#include<vector>
namespace ofec {
	template<typename TInterpreter>
	class DistanceCalculator  {
	public:

		using SolutionType = typename TInterpreter::SolutionType;
		using InterpreterType = typename TInterpreter;
	protected:

		SolutionType m_center;

		std::vector<std::vector<Real>> m_edge_weight;
		std::vector<std::vector<Real>> m_cur_edge_weight;
		std::vector<std::vector<int>> m_cur_edge_visited;
		Real m_curEdgeAvgDis = 0;
		Real m_curPosAvgDis = 0;
		Real m_curEdgeBestDis = 0;
		Real m_curPosBestDis = 0;
		// relatively ratio of the current edge weight
		Real m_alpha = 0.1;
	public:
		virtual void initialize(Algorithm *alg,Problem *pro) {
			for (auto& it : m_edge_weight) {
				std::fill(it.begin(), it.end(), 0);
			}
		}
		virtual void resize(const InterpreterType& interpreter) {
			auto& matSize(interpreter.getMatrixSize());
			m_edge_weight.resize(matSize.size());
			for (int idx(0); idx < m_edge_weight.size(); ++idx) {
				int from(m_edge_weight[idx].size());
				m_edge_weight[idx].resize(matSize[idx]);
				std::fill(m_edge_weight[idx].begin()+from, m_edge_weight[idx].end(), 0);
			}
			m_cur_edge_weight.resize(matSize.size());
			for (int idx(0); idx < m_edge_weight.size(); ++idx) {
				m_cur_edge_weight[idx].resize(matSize[idx]);
			}
			m_cur_edge_visited.resize(matSize.size());
			for (int idx(0); idx < m_edge_weight.size(); ++idx) {
				m_cur_edge_visited[idx].resize(matSize[idx]);
			}
		}
		void addDistance(Problem *pro,const InterpreterType& interpreter, 
			const std::vector<SolutionType*>& pop, 
			const std::vector<Real> & weight) {
			int centerIdx(0);
			for (int idx(0); idx < weight.size(); ++idx) {
				if (weight[idx] > weight[centerIdx]) {
					centerIdx = idx;
				}
			}
			for (auto& it : m_cur_edge_weight) {
				std::fill(it.begin(), it.end(),0);
			}
			for (auto& it : m_cur_edge_visited) {
				std::fill(it.begin(), it.end(), 0);
			}

			//std::cout << "popSize\t" << pop.size() << std::endl;
			for (int popId(0); popId < pop.size(); ++popId) {
			//	std::cout << "popId\t" <<  popId << std::endl;
				interpreter.updateEdges(pro, *pop[popId]);
				for (auto& edge : pop[popId]->edges()) {
					m_cur_edge_weight[edge.first][edge.second] += weight[popId];
					++m_cur_edge_visited[edge.first][edge.second];
				}
			//	std::cout << "popId\t" << popId << std::endl;
			}

			for (int x(0); x < m_cur_edge_visited.size(); ++x) {
				for (int y(0); y < m_cur_edge_visited[x].size(); ++y) {
					if (m_cur_edge_visited[x][y]) {
						m_cur_edge_weight[x][y] /= Real(m_cur_edge_visited[x][y]);
					}
					else m_cur_edge_weight[x][y] = 0;
					m_edge_weight[x][y] = m_edge_weight[x][y] * (1 - m_alpha) + m_cur_edge_weight[x][y]*m_alpha;
				}
			}
			m_curEdgeAvgDis = 0;
			for (auto& popIter : pop) {
				m_curEdgeAvgDis += edgeDistance(pro, interpreter, *popIter);
			}
			m_curEdgeAvgDis /= pop.size();
	        	
			m_center = *pop[centerIdx];
			m_curEdgeBestDis = edgeDistance(pro, interpreter, m_center);

			m_curPosAvgDis = 0;
			for (auto& popIter : pop) {
				m_curPosAvgDis += posDistance(pro, *popIter);
			}
			m_curPosAvgDis /= pop.size();

			m_curPosBestDis = posDistance(pro, m_center);
		}
		Real edgeDistance(Problem *pro, const InterpreterType& interpreter, SolutionType& cur) {
			interpreter.updateEdges(pro, cur);
			Real weight(0);
			for (auto& edge : cur.edges()) {
				weight += m_edge_weight[edge.first][edge.second];
			}
			return weight;
		}
		Real posDistance(Problem *pro,SolutionType& cur) {
			return m_center.variableDistance(cur, pro);
		}

		void printInfo() {
			std::cout << "avgEdgeDis\t" << m_curEdgeAvgDis << "\tBestEdgeDis\t" << m_curEdgeBestDis << "\t";
			std::cout << "avgPosDis\t" << m_curPosAvgDis << "\tBestPosDis\t" << m_curPosBestDis << std::endl;
		}

	};
}


#endif