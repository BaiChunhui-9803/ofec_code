#ifndef OFEC_NBN_EAX_POP_V4_V2_H
#define OFEC_NBN_EAX_POP_V4_V2_H

#include "../eax_nbn_4/nbn_eax_pop_v4.h"

#include "../../../../../utility/hnsw/hnsw_solution/hnsw_sol_model.h"
#include "../../../../../utility/nbn_visualization/nbn_nearest_better_calculator.h"


#include "../../../../../utility/nbn_visualization/tree_graph_simple.h"

namespace ofec {




	class PopNBN_EAX_V4_V2 : public PopNBN_EAX_V4 {
	protected:
		double m_searchRadius = 50;
		int m_maxStagnationTime = 30;
		int neighborSize = 50;
		bool m_expandStage = true;
		n2::HnswSolModel m_hnswModel;


		//std::vector<int> m_dividison_basin;
		//std::vector<int> m_search_basin;

		std::vector<unsigned> m_solCreateIter;
		std::vector<unsigned> m_solUpdateIter;

		std::vector<double> m_dis2parent;

		Indi* m_curPeak = nullptr;

		//std::vector<int> m_belong;

		//std::map<unsigned long long, eax_tsp::TIndi*>  m_solMap;
		//NBN_hash m_solHash;
		int evolveWithHnsw(Problem* pro, Algorithm* alg, Random* rnd);

		//bool judgeSolByClass(Indi* cur, int curClass  );
		void generateSolRandom(ofec::Problem* pro);

		void calNodeState(Indi* cur);
		
		void getParentTwoNeighbor(Indi* cur, SolutionBase*& parent, ofec::Problem* pro);



		bool judgePeak(Indi* cur, double dis2parent) {

			if (cur->m_explorationTimes == 0 && dis2parent > m_searchRadius / 2.0) {
				if (cur->m_state == Indi::NodeState::kUnknow) {
					calNodeState(cur);
				}
				return cur->m_state == Indi::NodeState::kPossiblePeak;
			}
			return false;
		//	return cur->m_explorationTimes==0&& dis2parent> m_searchArea/2.0&&
		}


		bool addSolToHistory(std::shared_ptr<Indi>& cursol, Indi*& indi) {
			auto& cursolx = cursol->variable().vect();
			indi = nullptr;
			auto hash = m_solHash.calHash(cursolx);
			indi = m_solMap[hash];
			if (indi == nullptr) {
				m_indis.emplace_back(std::move(cursol));
				m_indis.back()->setRndId(hash);
				//	m_indis.back()->setCurClass(0);
				indi = m_indis.back().get();
				//		m_curPop.push_back(indi);
				m_indis.back()->setSolId(m_hnswModel.AddData(indi));
				return true;
			}
			else return false;
			//return indi;
		}
	public:

		void calNBN(
			std::vector<TreeGraphSimple::NBNdata>& nbn_data,
			//std::vector<double>& dataBasin,
			std::vector<TreeGraphSimple::NBNdata>& nbn_data_best,
			//std::vector<double> 
			std::vector<int>& curPop,
			std::vector<int>& curBasin,
			std::vector<int>& colorfulArea,
			std::vector<int>& bestIds, 
			ofec::Problem* pro, 
			ofec::Random* rnd
			);

		virtual void initialize(Problem* pro, Random* rnd) override;
	//	void expandSearch();

		int evolve(Problem* pro, Algorithm* alg, Random* rnd) override {
			return evolveWithHnsw(pro, alg, rnd);
		}

		
	};
}



#endif