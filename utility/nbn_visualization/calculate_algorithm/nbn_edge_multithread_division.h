#ifndef NBN_EDGE_MULTITHREAD_DIVISION_H
#define NBN_EDGE_MULTITHREAD_DIVISION_H

#include <vector>
#include <memory>
#include "../../../core/problem/solution.h"
#include "../../../core/random/newran.h"
#include "../../../utility/hash_table/hash_table.h"
#include "nbn_division_base.h"


//#include "core/problem/solution.h"
//#include "utility/typevar/typevar.h"

// for edge tsp 

namespace ofec {
	class NBN_EdgeMultiThreadDivision : public NBN_DivisionBase {
		struct NetworkNode {
			//int m_network_id = 0;
			std::array<int, 2> m_network_id = { -1,-1 };
			std::vector<NetworkNode*> m_parents;
		//	std::vector<NetworkNode*> m_sons;
			std::vector<int> m_repre_edge;
		};

		struct NetworkNodeSolBatch {
			std::vector<int> m_indis_ids;
			int m_best_solId = -1;
		};

		std::vector<std::vector<NetworkNode>> m_network;
		//std::vector<std::vector<NetworkNodeSolBatch>> networkSols;



		// parameters 
		int m_divide_threadhold = 2e2;
		int m_edge_division = 50;
		int m_num_division = 1e3;
		int m_divide_times = 0;
//		int m_dim = 0;
		int m_num_loop = 5;

		int m_numSubIndis = 5;

		int m_numSubLoop = 10;
		double m_err = 0.1;
		//std::vector<int> m_edge_seq;

		//std::vector<int> 
		//std::vector<int> m_selected_edges;
		//std::vector<std::vector<NetworkNode>> m_networks;
		////std::vector<std::vector<int>> ;
		//std::vector<std::vector<int>> m_edgeToEdgeId;
		//std::vector<int> m_edge_indis_ids;
		//std::vector<int> m_sol_direction;
		std::vector<std::vector<std::array<int, 2>>> m_sols_edges;
		std::vector<std::shared_ptr<SolutionBase>> m_sols;
		std::vector<double> m_fitness;
		std::vector<int> m_belong;
		std::vector<double> m_dis2parent;
		std::vector<bool> m_flagOpt;

		int m_bestSolId = -1;

		std::vector<int> m_sortedIds;
		

		int m_centerSolId = -1;
		std::vector<std::array<int, 2>> m_centerSolEdge;


	protected:


		bool compareSol(int a, int b)const {
			if (m_fitness[a] == m_fitness[b]) {
				for (int idx(0); idx < m_sols_edges[a].size(); ++idx){
					for (int idy(0); idy < 2; ++idy) {
						if (m_sols_edges[a][idx][idy] != m_sols_edges[b][idx][idy]) {
							return m_sols_edges[a][idx][idy] < m_sols_edges[b][idx][idy];
						}
					}
				}

				return a < b;
			}
			else return m_fitness[a] > m_fitness[b];
		}


		void resize(int size) {
			//m_sol_direction.resize(size);
			m_sols_edges.resize(size);
			m_sols.resize(size);
			m_fitness.resize(size, 0);
			m_belong.resize(size, -1);
			m_dis2parent.resize(size, 1.0);
			m_flagOpt.resize(size, false);
		}

		void swapSol(int idx, int idy) {
			std::swap(m_sols_edges[idx], m_sols_edges[idy]);
			std::swap(m_sols[idx], m_sols[idy]);
			std::swap(m_fitness[idx], m_fitness[idy]);
			//swap(m_fitness[idx], m_fitness[idy]);
			std::swap(m_belong[idx], m_belong[idy]);
			std::swap(m_dis2parent[idx], m_dis2parent[idy]);

			bool a = m_flagOpt[idx];
			m_flagOpt[idx] = m_flagOpt[idy];
			m_flagOpt[idy] = a;
		//	std::swap<bool>(m_flagOpt[idx], m_flagOpt[idy]);
		}

		void solToIdx(int solId,
			const std::vector<int> & edge_seq, 
			const std::vector<int> & sol_direction,
			int& nodeId)const;
		void idxToVec(int idx, std::vector<int>& cur)const;
		void vecToIdx(const std::vector<int>& cur, int& idx)const;

		void Cnm(std::vector<std::vector<bool>>& selectedIdx,
			std::vector<bool>& cur, int m, int from)const;
		void Cnm(std::vector<unsigned>& selectedIdx, unsigned cur, int curSize, int m, int from)const;


		



		void generateSolsThreadTask(int startId, int from, int to, double seed);
		void generateSol(int solId, Random *rnd);

		void updateSolsThreadTask(int startedId, const std::vector<std::shared_ptr<SolutionBase>>& sols,int from, int to);
		void updateSol(int solId, const std::vector<std::shared_ptr<SolutionBase>>& sols);


		void calHashThreadTask(int from, int to,
			std::vector<unsigned>& hashVal, const utility::HashRandomValue& hash_table);
	



	/*	void filterSolThreadTask(int from, int to,
			std::vector<bool>& flagSame, 
		);*/
		//	void updateNBN()const;
		// update solutions neighbor
		void clearNetworkSols(std::vector<std::vector<NetworkNodeSolBatch>>& network_sols)const;
		void divideByEdge(NetworkNodeSolBatch& node,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			Random *rnd,
			const std::vector<int> & edge_seq,
			const std::vector<int>& sol_dir,
			int divDim, bool divFlag)const;
		// return the best idxs
		int udpateNeighbor(std::vector<int>& indis,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			Random *rnd
			)const;


		void updateInfo(int curId,int betterId,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			Random* rnd
		)const;
		
		void updateDivisionTask(
			double seed,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			int numLoop = 5
		)const;


		void updateDivisionSubRegionTask(
			double seed,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			int numLoop = 5
		)const;



		void updateDivisionSubRegionBandEdgeTask(
			double seed,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			int numLoop = 5
		)const;


		int updateDivisionTaskSubRegionDivisionByEdge(
			ofec::Random* rnd,
			std::vector<int>& selectedEdgeIds, 
			int from,
			const std::vector<int>& solIds,
			std::vector<int>& belong,
			std::vector<double>& dis2parent
		)const;


		int updateDivisionTaskSubRegionDivisionByEdgeBandEdge(
			ofec::Random* rnd,
			std::vector<int>& selectedEdgeIds,
			int from,
			const std::vector<int>& solIds,
			std::vector<int>& belong,
			std::vector<double>& dis2parent
		)const;
		


		int updateDivisionTaskSubRegionDivisionByEdgeOnCenterEdge(
			ofec::Random* rnd,
			std::vector<int>& selectedEdgeIds,
			int from,
			const std::vector<int>& sortedSolIds,
			std::vector<int>& belong,
			std::vector<double>& dis2parent
		)const;


		void updateDivisionTaskSubRegionDivisionByEdgeOnCenterEdgeTask(
			double seed,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			int numLoop = 5
		)const;


		double solDis(int solIdx, int solIdy)const;

		void mergeBelongInfoThreadTask(
			double seed, int from, int to,
			const std::vector<std::vector<int>>& total_belong,
		    const std::vector<std::vector<double>> total_dis2parent
/*,		    std::vector<int>& belong,
			std::vector<double>& dis2parent*/);

		

	public:

		void getResult(std::vector<int>& belong,std::vector<double>& dis2parent) {
			belong = m_belong;
			dis2parent = m_dis2parent;
		}


		void updateBestSolSub(int numSubSols) {

			//m_bestSolId = numSubSols;
			
			std::vector<int> solIds(m_sols.size());
			for (int idx(0); idx < solIds.size(); ++idx) {
				solIds[idx] = idx;
			}
			std::sort(solIds.begin(), solIds.end(), [&](int a, int b) {
				return compareSol(a, b);
				//	return m_fitness[a] > m_fitness[b];
			});

			solIds.resize(numSubSols);
			m_bestSolId = udpateNeighbor(solIds, m_belong, m_dis2parent, m_random.get());

		}


		void updateSolsDis2Sol(const std::vector<int>& division, int bestSolId,
			std::vector<double>& dis2parent, std::vector<int>& belong, ofec::Random* rnd
			) const {
			for (auto& solId : division) {
				if (solId != bestSolId) {
					double curDis = solDis(solId, bestSolId);
					if (curDis < dis2parent[solId]) {
						dis2parent[solId] = curDis;
						belong[solId] = bestSolId;
					}
					else if (curDis == m_dis2parent[solId] && rnd->uniform.next() < 0.5) {
						belong[solId] = bestSolId;
					}
				}
			}
		}

		void updateSolsDis2Center() {
			for (int solId(0); solId < m_sols.size(); ++solId) {
				if (solId != m_centerSolId) {
					if(compareSol(m_centerSolId, solId)) {
						double curDis = solDis(solId, m_centerSolId);
						if (curDis < m_dis2parent[solId]) {
							m_dis2parent[solId] = curDis;
							m_belong[solId] = m_centerSolId;
						}
						else if (curDis == m_dis2parent[solId] && m_random->uniform.next() < 0.5) {
							m_belong[solId] = m_centerSolId;
						}
					}
				}
			}
		}


		void setCenterSol(int centerSolId);

		void setCenterSol(const std::shared_ptr<ofec::SolutionBase>& centerSol);


		void updateNetwork();

		void setNumberLoop(int numLoop) {
			m_num_loop = numLoop;
		}
		void updateDivision();
		void updateDivisionMultithread();
		void updateDivisionSubRegionMultiThread();

		void updateDivisionSubRegionBandEdgeMultiThread();
		void updateDivisionSubRegionOnCenteredEdgeMultiThread();

		void updateSols(const std::vector<std::shared_ptr<SolutionBase>>& sols);
		void generateSols(int numSample);
		
		void filterSameSolutions(int from,int to);

		// Í¨¹ý NBN_DivisionBase ¼Ì³Ð
		virtual void initialize_(bool flag_grid_sample = true) override;
		virtual size_t size() const override;
		virtual void addSol(const SolutionBase& new_sol, int& belong_id, bool flag_opt = false, int popIter = -1, int popSolId = -1, int algId = -1) override;
		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<bool>& flagOpt) const override;
		virtual void getSharedNBN(std::vector<std::shared_ptr<SolutionBase>>& sols, std::vector<double>& fitness, std::vector<int>& belong, std::vector<double>& dis2parent, std::vector<int> popIters, std::vector<int> popSolIds, std::vector<int> algIds) const override;

	};
}

#endif