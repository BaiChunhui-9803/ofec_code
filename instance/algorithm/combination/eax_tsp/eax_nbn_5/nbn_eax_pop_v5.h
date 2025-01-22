#ifndef OFEC_NBN_EAX_POP_V5_H
#define OFEC_NBN_EAX_POP_V5_H

#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/problem/solution.h"

#include "../../../../algorithm/realworld/DVRP/LKH/INCLUDE/LKH.h"

#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../environment.h"
#include "../../nbn_alg_com/gl_calculator.h"
#include <queue>
#include <set>

namespace ofec {

	class PopNBN_EAX_V5 : public Population<Solution<VarVec<int>>>, public eax_tsp::TEnvironment {
	//public:

		struct NBN_hash {
			std::vector<unsigned long long> m_hash_val;
			void initialize(ofec::Random* rnd, int numVariable) {
				m_hash_val.resize(numVariable);
				unsigned long long maxVal = std::numeric_limits<unsigned long long>::max();
				maxVal = sqrt(maxVal);
				for (auto& it : m_hash_val) {
					it = rnd->uniform.nextNonStd<unsigned long long>(0, maxVal);
				}

			}
			unsigned long long getHash(int idx) {
				return m_hash_val[idx];
			}
			unsigned long long calHash(const std::vector<int>& sol)const {
				unsigned long long Hash(0);
				for (int idx(1); idx < sol.size(); ++idx) {
					Hash ^= m_hash_val[sol[idx - 1]] * m_hash_val[sol[idx]];
				}
				Hash ^= m_hash_val[sol.back()] * m_hash_val[sol.front()];
				return Hash;
			}
		};




		struct  RelationList;
		class Indi;


		struct RelationShip {
			Indi* m_cur = nullptr;
			Indi* m_other = nullptr;
			double m_distance = std::numeric_limits<double>::max();
			
			void initialize() {
				m_cur = nullptr;
				m_other = nullptr;
				m_distance = std::numeric_limits<double>::max();
			}
			void set(Indi* cur, Indi* parent, double dis = 0) {
				m_cur = cur;
				m_other = parent;
				m_distance = dis;
			}

			inline bool isActive()const {
				return m_cur != nullptr;
			}
		};

		struct RelationNode {
			RelationShip m_cur;
			RelationNode* m_before = nullptr;
			RelationNode* m_after = nullptr;
			RelationList* m_list = nullptr;

			
			inline double getDis() const{
				return m_cur.m_distance;
			}
			void initialize() {
				m_cur.initialize();
				m_before = nullptr;
				m_after = nullptr;
				m_list = nullptr;
			}
			inline bool isActive()const {
				return m_cur.isActive();
			}

			void set(Indi* cur, Indi* parent, double dis = 0) {
				m_cur.set(cur, parent, dis);
			}


			inline Indi* getCur();
			inline Indi* getParent();


			void removeFromLink();
			void insertToLink(RelationList* list);
		};

		struct RelationList {
			std::shared_ptr<RelationNode> m_head;
			std::shared_ptr<RelationNode> m_tail;
			int m_num = 0;
			//double m_sumDis2parent = 0;


			RelationList() {
				m_head.reset(new RelationNode);
				m_tail.reset(new RelationNode);
			}

			void initialize(Indi* cur) {
				m_head->initialize();
				m_head->set(cur, nullptr);
				m_tail->initialize();
				m_tail->set(cur, nullptr);
				m_head->m_after = m_tail.get();
				m_tail->m_before = m_head.get();
				m_num = 0;
			}

			void insertNode(RelationNode* cur);
			void removeNode(RelationNode* cur);

			void clearLink();
		};



		class Indi : public eax_tsp::TIndi {
		protected:

		public:

			int m_solId = -1;

			unsigned long long m_randId = 0;

			RelationNode* m_direct_parent;
			std::vector<std::shared_ptr<RelationNode>> m_parents;
			RelationList m_sons;
		//	double m_fitness = 0;


			unsigned m_createdIter = 0;
			int m_stagnation_time = 0;
			bool m_improve = false;
			unsigned long long m_curClassFlag = 0;
			unsigned long long m_curVisitedId = 0;

			Indi() {}
			~Indi() = default;

			void updateState() {
				m_stagnation_time = 0;
				m_improve = false;
			}

			void initialize(int N, int solId, unsigned createdIter) {
				define(N);
				m_solId = solId;
				m_randId = 0;
				m_direct_parent = nullptr;
				m_parents.clear();
				m_sons.initialize(this);
				m_fitness = 0;
				m_createdIter = 0;
				m_stagnation_time = 0;
				m_improve = false;
				m_curClassFlag = 0;
				m_curVisitedId = 0;
			}

			void setSolId(int solId) {
				m_solId = solId;
			}
			void setRndId(unsigned long long rndId) {
				m_randId = rndId;
			}
			void setCurClass(unsigned long long curClass) {
				m_curClassFlag = curClass;
			}
			void copySol(const Indi& other) {
				fLink = other.fLink;
				fEvaluationValue = other.fEvaluationValue;
			}

			void clearParent() {
				if (m_direct_parent != nullptr) m_direct_parent->removeFromLink();
				m_parents.clear();
				m_direct_parent = nullptr;
			}




			bool updateParent(Indi* parent, double minDis, ofec::Random* rnd) {
				bool updateFlag = false;
				if (m_parents.empty()) {
					updateFlag = true;
				}
				else {
					if (m_parents.front()->getDis() == minDis) {
						updateFlag = true;
					}
					else if (m_parents.front()->getDis() > minDis) {
						m_parents.clear();
						updateFlag = true;
					}
				}

				if (updateFlag) {
					if (m_direct_parent != nullptr) m_direct_parent->removeFromLink();
					m_parents.emplace_back(new RelationNode());
					m_parents.back()->initialize();
					m_parents.back()->set(this, parent, minDis);
					m_direct_parent = m_parents[rnd->uniform.nextNonStd<int>(0, m_parents.size())].get();
					parent->m_sons.insertNode(m_direct_parent);

				}

			}

			void setRoot(Indi* root) {

				double minDis = std::numeric_limits<double>::max();
				m_parents.emplace_back(new RelationNode());
				m_parents.back()->initialize();
				m_parents.back()->set(this, root, minDis);
				m_direct_parent = m_parents.front().get();
				root->m_sons.insertNode(m_direct_parent);

			}

		};


	public:
		PopNBN_EAX_V5() : Population(), TEnvironment() {}
		PopNBN_EAX_V5(size_t size_pop, Problem* pro) : Population(size_pop, pro), TEnvironment() {};

		virtual ~PopNBN_EAX_V5() = default;
	
	protected:

		unsigned m_curIter = 0;
		double m_pos = 1.0;
		unsigned m_numPop = 100;
		unsigned m_maxStag = 50;
		unsigned m_maxSampleTime = 1e3;
		const unsigned long long m_maxULL = std::numeric_limits<unsigned long long>::max();
		unsigned long long m_maxCurClassFlag = 0;
		unsigned long long m_curVisitedTime = 0;

		Indi* m_peakAround = nullptr;

		unsigned long long m_curEvolveClass = 0;
		std::vector<Indi*> m_curPop;
		std::vector<Indi*> m_stagnationIndis;
		std::vector<std::shared_ptr<Indi>> m_indis;
		std::shared_ptr<Indi> m_root;

		std::vector<int> mc_StagnationTimes;
		int m_curMaxStag = 0;


		std::vector<Indi*> m_diverSols;




		Indi* m_curPeak;



		std::vector<Indi*> m_nondiminatedStates;


		std::map<unsigned long long, Indi*>  m_solMap;
		NBN_hash m_solHash;


		//	std::vector<int> m_id2curId;
		std::vector<Indi*> m_peaks;
		std::vector<eax_tsp::TIndi> m_bestSols;
		double m_bestObj = 0;
		bool m_expanedStage = false;
		LKH::LKHAlg m_lkh_alg;
		int m_moveType = 5;

		void updateFitness(Indi* indi) {
			indi->setFitness( m_pos * indi->fEvaluationValue);
		}

		//void getSons();



	public:

		void getPopIds(std::vector<int>& popIds);


		void calNBN(
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<double>& fitness,
			std::vector<int>& optNBNid,
			std::vector<double>& optFit,
			ofec::Random* rnd
		);

	};
}

#endif // !OFEC_PopGL_NBN_COM_ALG_H