#ifndef OFEC_NBN_EAX_POP_V4_H
#define OFEC_NBN_EAX_POP_V4_H

#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/problem/solution.h"

#include "../../../../algorithm/realworld/DVRP/LKH/INCLUDE/LKH.h"

#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../environment.h"
#include "../../nbn_alg_com/gl_calculator.h"
#include <queue>
#include <set>
namespace ofec {


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



		//unsigned long long calHash(
		//	const std::vector<std::vector<int>>& solLink)const {
		//	unsigned long long Hash(0);
		//	for (int idx(0); idx < m_hash_val.size(); ++idx) {
		//		Hash ^= m_hash_val[idx] * m_hash_val[solLink[idx].front()];
		//	}
		//	return Hash;
		//}

		

	};


	class PopNBN_EAX_V4 : public Population<Solution<VariableVector<int>>>, public eax_tsp::TEnvironment {


	protected:


		struct  RelationList;
		class Indi;
		
		struct RelationNode {
			Indi* m_cur;
			RelationNode* m_before = nullptr;
			RelationNode* m_after = nullptr;
			RelationList* m_list = nullptr;

			
			void initialize(Indi* cur) {
				m_cur = cur;
				m_before = nullptr;
				m_after = nullptr;
				m_list = nullptr;
			}

			void set(Indi* cur) {
			    m_cur=cur;
			}

			inline Indi* getCur();
			inline Indi* getParent();
			

			void removeFromLink();
			void insertToLink(RelationList* list);
		};

		struct RelationList {
			std::shared_ptr<RelationNode> m_head;
			std::shared_ptr<RelationNode> m_tail;
			Indi* m_cur = nullptr;

			int m_num = 0;
			//double m_sumDis2parent = 0;


			RelationList() {
				m_head.reset(new RelationNode);
				m_tail.reset(new RelationNode);
			}

			void initialize(Indi* cur) {
				m_head->initialize(nullptr);
				m_tail->initialize(nullptr);
				m_head->m_after = m_tail.get();
				m_tail->m_before = m_head.get();
				m_cur = cur;
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
			int m_explorationTimes = 0;
			unsigned long long m_curClassFlag = 0;
			unsigned long long m_searchArea = 0;
			enum class  NodeState
			{
				kUnknow, kNormal, kPossiblePeak
			};
			NodeState m_state = NodeState::kUnknow;
			





			bool m_exploitationFlag = false;
			

			unsigned long long m_randId = 0;
			std::shared_ptr<RelationNode> m_parentNode;
			RelationList m_sons;
			

			double m_originDis = 0;
		
			double m_dis2parent = std::numeric_limits<double>::max();
		//	double m_fitness = 0;
			
			unsigned m_createdIter = 0;
			int m_stagnation_time = 0;
			bool m_improve = false;
			bool m_stag = false;

			unsigned long long m_curVisitedId = 0;
			
			bool m_flagMultiParent = false;
			

			bool m_flagOpt = false;
			

			double m_maxRadius = 0;

			double m_curSearchDis = 0;
			

			Indi() {
				m_parentNode.reset(new RelationNode);
			}

			~Indi() = default;

			void updateState() {
				m_stagnation_time = 0;
				m_improve = false;
				m_stag = false;
			}

			void initialize(int N, int solId, unsigned createdIter) {
				define(N);
				m_solId = solId;
				m_createdIter = createdIter;

				m_parentNode->initialize(this);
				m_dis2parent = std::numeric_limits<double>::max();
				m_fitness = 0;
				m_stagnation_time = 0;
				m_curClassFlag = 0;
				m_sons.initialize(this);
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


			void clearParent() {
				m_parentNode->removeFromLink();
	//			m_parentNode->m_cur = nullptr;
			}


			void setMultiParent(bool flag) {
				m_flagMultiParent = flag;
			}

			void updateParent(Indi* parent, double minDis) {
				parent->m_sons.insertNode(m_parentNode.get());
				if (m_dis2parent > minDis) {
					m_dis2parent = minDis;
					m_stagnation_time = 0;
				}
				//if (m_dis2parent<0) {
				//	int stop = -1;
				//}
			}
			void updateParent(Indi* parent) {
		
				parent->m_sons.insertNode(m_parentNode.get());
				//m_dis2parent = minDis;
				//m_stagnation_time = 0;
			}

			void setRoot(Indi* root) {
				root->m_sons.insertNode(m_parentNode.get());

			}


			int judgeInside(Indi* parent, double dis) {
				if (parent->m_fitness > m_fitness) {
					if (dis == m_dis2parent) {
						return 0;
					}
					else if (dis < m_dis2parent) {
						return 1;
					}
				}
				return -1;
			}
			// a<b: -2, -1, 0, a>b 1, 2
			static int judgeRelationship(
				const Indi* a,
				const Indi* b,
				double& dis
			) {
				dis = -1;
				if (a->m_fitness < b->m_fitness) {
					dis = a->distanceTo(*b);
					if (dis == a->m_dis2parent) {
						return -1;
					}
					else if (dis < a->m_dis2parent) {
						return -2;
					}
				}
				else if (a->m_fitness > b->m_fitness) {
					dis = a->distanceTo(*b);
					if (dis == b->m_dis2parent) {
						return 1;
					}
					else if (dis < b->m_dis2parent) {
						return 2;
					}
				}


				else if (a->m_solId < b->m_solId) {
					dis = a->distanceTo(*b);
					if (dis == a->m_dis2parent) {
						return -1;
					}
					else if (dis < a->m_dis2parent) {
						return -2;
					}
				}
				else {
					dis = a->distanceTo(*b);
					if (dis == b->m_dis2parent) {
						return 1;
					}
					else if (dis < b->m_dis2parent) {
						return 2;
					}
				}
				return 0;
			}

			static int judgeRelationshipWithDis(
				const Indi* a,
				const Indi* b,
				double dis
			) {
				if (a->m_fitness < b->m_fitness) {
					if (dis == a->m_dis2parent) {
						return -1;
					}
					else if (dis < a->m_dis2parent) {
						return -2;
					}
				}
				else if (a->m_fitness > b->m_fitness) {
					if (dis == b->m_dis2parent) {
						return 1;
					}
					else if (dis < b->m_dis2parent) {
						return 2;
					}
				}
				else if (a->m_solId < b->m_solId) {
					if (dis == a->m_dis2parent) {
						return -1;
					}
					else if (dis < a->m_dis2parent) {
						return -2;
					}
				}
				else {
					if (dis == b->m_dis2parent) {
						return 1;
					}
					else if (dis < b->m_dis2parent) {
						return 2;
					}
				}
				return 0;
			}


			static int judgeRelationshipFit(
				const Indi* better,
				const Indi* worse,
				double dis
			) {
				{
					if (dis == worse->m_dis2parent) {
						return 1;
					}
					else if (dis < worse->m_dis2parent) {
						return 2;
					}
				}
				return 0;
			}



			static void updateRelationship(
				Indi* a,
				Indi* b,
				int rf, double curdis,
				ofec::Random * rnd) {
				if (rf == 1) {
					b->setMultiParent(true);
					if (rnd->uniform.next() < 0.5) {
						b->updateParent(a, curdis);
					}
				}
				else if (rf == 2) {
					b->updateParent(a, curdis);
				}
				else if (rf == -1) {
					a->setMultiParent(true);
					if (rnd->uniform.next() < 0.5) {
						a->updateParent(b, curdis);
					}
				}
				else if (rf == -2) {
					a->updateParent(b, curdis);
				}
			}




		};


		bool generateSols(Indi* cur, 
			std::shared_ptr<Indi>& newSol, 
			ofec::Problem* pro, 
			ofec::Random* rnd);
		

		void generateRandomSols(Indi* cur, 
			std::vector<Indi*>& newsols, 
			int numSamples, int maxSamples,
			ofec::Problem* pro, ofec::Random* rnd);
		

		void anayalizeNewSols(std::vector<Indi*>& newsols, ofec::Random* rnd);

	public:
		PopNBN_EAX_V4() : Population(), TEnvironment() {}
		PopNBN_EAX_V4(size_t size_pop, Problem* pro) : Population(size_pop, pro), TEnvironment() {};

		virtual ~PopNBN_EAX_V4() = default;
		int evolve2(Problem* pro, Algorithm* alg, Random* rnd);
		int evolve3(Problem* pro, Algorithm* alg, Random* rnd);

		int evolve(Problem* pro, Algorithm* alg, Random* rnd) override;


		
		
		virtual void initialize(Problem* pro, Random* rnd) override;

		void increaseClassFlag() {
			if (++m_maxCurClassFlag == m_maxULL) {
				for (auto& it : m_indis) {
					it->m_curClassFlag = 0;
				}
				m_root->m_curClassFlag = 0;
				m_maxCurClassFlag = 1;
			}
		
		}
		void increaseVisitedTimes() {
			if (++m_curVisitedTime == m_maxULL) {
				for (auto& it : m_indis) {
					it->m_curVisitedId = 0;
				}
				m_root->m_curVisitedId = 0;
				m_curVisitedTime = 1;
			}
		}


		void generateSols(
		std::vector<std::shared_ptr<Indi>>& newSols,
		ofec::Problem* pro, ofec::Random* rnd) {
			std::shared_ptr<Indi> curindi;
			std::vector<int> cursol;
			int dim = pro->numberVariables();
			while (m_curPop.size() + newSols.size() < m_numPop) {
				curindi.reset(new Indi);
				curindi->initialize(dim, -1, m_curIter);
				auto& randIndi = m_curPop[rnd->uniform.nextNonStd<int>(0, m_curPop.size())];

				randIndi->transferSol(cursol);
			//	int solIdx = m_id2curId[randIndi->m_solId];
				int radius = std::round(rnd->uniform.nextNonStd<int>(1.0, randIndi->m_maxRadius / 4.0));
				if (radius) {

					while (radius--) {
						auto a = rnd->uniform.nextNonStd<int>(0, dim);
						auto b = rnd->uniform.nextNonStd<int>(0, dim);
						std::swap(cursol[a], cursol[b]);
					}
					newSols.push_back(curindi);


				}
			}


		}

		void getBorder(const std::vector<Indi*>& curIndi,
			std::vector<Indi*>& border) {
			increaseVisitedTimes();
			for (auto& it : curIndi) {
				it->m_curVisitedId = m_curVisitedTime;
				//	neighbors.push_back(it);
			}


			for (auto& it : curIndi) {
				auto parent = it->m_parentNode->getParent();
				if (parent->m_solId != -1){
					if (parent->m_curVisitedId != m_curVisitedTime) {
						parent->m_curVisitedId = m_curVisitedTime;
						border.push_back(parent);
					}
				}
			
				RelationNode* head = parent->m_sons.m_head.get();
				while (head->m_after->m_cur != nullptr) {
					head = head->m_after;
					if (head->m_cur->m_curVisitedId != m_curVisitedTime) {
						head->m_cur->m_curVisitedId = m_curVisitedTime;
						border.push_back(head->m_cur);
					}
				}

			}
			
		}

		void getNeighbors(
			const std::vector<Indi*>& curIndi,
			std::vector<Indi*>& neighbors
		) {
			increaseVisitedTimes();
			std::queue<Indi*> que;
			for (auto& it : curIndi) {
				it->m_curVisitedId = m_curVisitedTime;
			//	neighbors.push_back(it);
			}

			for (auto& it : curIndi) {
				auto parent = it->m_parentNode->getParent();
				if (parent!=nullptr&&parent->m_curVisitedId != m_curVisitedTime) {
					que.push(it);
					parent->m_curVisitedId = m_curVisitedTime;
					if (parent->m_solId != -1)
						neighbors.push_back(it);
				}
				
			}

			
			while (!que.empty()) {
				auto it = que.front();
				que.pop();
				{
					RelationNode* head = it->m_sons.m_head.get();
					while (head->m_after->m_cur != nullptr) {
						head = head->m_after;
						if (head->m_cur->m_curVisitedId != m_curVisitedTime) {
							neighbors.push_back(head->m_cur);
							head->m_cur->m_curVisitedId = m_curVisitedTime;

						}
					}
				}

				auto parent = it->m_parentNode->getParent();
				if (parent != nullptr && parent->m_curVisitedId != m_curVisitedTime) {
					que.push(parent);
					parent->m_curVisitedId = m_curVisitedTime;

					if (parent->m_solId != -1)
						neighbors.push_back(parent);			
				}
			}
		}

		bool judgeInside(std::queue<Indi*>& que, 
			Indi* curindi, bool flagInside, unsigned long long curClass);
		

		void updateIndiActiveFlag(
			std::vector<bool>& actives,
			std::vector<Indi*>& newIndi
		) {
			for (int idx(0); idx < newIndi.size(); ++idx) {
				if (newIndi[idx]->m_dis2parent == 0) {
					actives[idx] = false;
				}
			}

		}

		void updateRoot(std::vector<Indi*>& newIndi){
			for (auto& it : newIndi) {

				if (it->m_parentNode->m_cur == nullptr) {
					it->setRoot(m_root.get());
				}
			}
		}

		void updateRelationShipOfV(
			std::vector<Indi*>& newIndi, ofec::Random* rnd) {
			double curdis(0);
			int rf(0);
			for (int idx(0); idx < newIndi.size(); ++idx) {
				auto a = newIndi[idx];
				for (int idy(idx + 1); idy < newIndi.size(); ++idy) {
					auto b = newIndi[idy];
					if (a->m_solId != b->m_solId) {
						rf = Indi::judgeRelationship(a, b, curdis);
						Indi::updateRelationship(a, b, rf, curdis, rnd);
					}
				}
			}

		//	updateIndiActiveFlag(actives, newIndi);
		}


		void updateRelationShipOf2V(
			std::vector<Indi*>& newIndi,
			std::vector<Indi*>& oldIndi,
			ofec::Random* rnd) {
			double curdis(0);
			int rf(0);
			//actives.resize(newIndi.size());
			for (int idx(0); idx < newIndi.size(); ++idx) {
				auto a = newIndi[idx];
				for (int idy(0); idy < oldIndi.size(); ++idy) {
					auto b = oldIndi[idy];
					if (a->m_solId != b->m_solId) {
						rf = Indi::judgeRelationship(a, b, curdis);
						Indi::updateRelationship(a, b, rf, curdis, rnd);
					}
				}
			}
			//updateIndiActiveFlag(actives, newIndi);
		}


		void insertBestSol(const std::string& filepath, ofec::Problem* pro);
		void updateNondiminatedStates(ofec::Problem* pro, ofec::Random* rnd);

		int exploitationState_v1(ofec::Problem* pro, ofec::Random* rnd);
		void getSubTree(Indi* cur, std::vector<Indi*>& subtrees);
		
		void getTotalSubTree(Indi* cur, std::vector<Indi*>& subtrees)const;
		
		void getDiversitySolInBasin(Indi* cur, std::vector<Indi*>& sols)const;

		void analysis(std::vector<int>& optNearIds);
		

		void setBasinsCurClass(Indi* cur, int curClass);


		Indi* addSolToIndis(std::shared_ptr<Indi>& cur,
			unsigned long long hashValue);


		void getDiverSols(std::vector<Indi*>& sols);
		

		void expandSearch_v1(Problem* pro, Algorithm* alg, Random* rnd);
		

		void expandSearch_v2(Problem* pro, Algorithm* alg, Random* rnd);
		
		void expandSearch_test(Problem* pro, Algorithm* alg, Random* rnd);
		void expandSearch_test2(Problem* pro, Algorithm* alg, Random* rnd);


		void filterNondominatedState(std::vector<Indi*>& totalStates)const;

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
		

		std::map<unsigned long long, Indi*>  m_solMap ;
		NBN_hash m_solHash;


	//	std::vector<int> m_id2curId;
		std::vector<Indi*> m_peaks;
		std::vector<eax_tsp::TIndi> m_bestSols;
		double m_bestObj = 0;
		bool m_expanedStage = false;
		LKH::LKHAlg m_lkh_alg;
		int m_moveType = 5;

		void updateFitness(eax_tsp::TIndi* indi) {
			indi->setFitness(m_pos * indi->fEvaluationValue);
		}

	//	void getSons();
		


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