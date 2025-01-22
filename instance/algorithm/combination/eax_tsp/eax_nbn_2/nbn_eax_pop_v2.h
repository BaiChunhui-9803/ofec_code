#ifndef OFEC_NBN_EAX_POP_V2_H
#define OFEC_NBN_EAX_POP_V2_H

#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/algorithm/individual.h"

#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../environment.h"
#include "../../nbn_alg_com/gl_calculator.h"

namespace ofec {




	class PopNBN_EAX_V2 : public Population<Individual<VarVec<int>>>, public eax_tsp::TEnvironment {
	
	public:
		class NBN_indi;
		struct  RelationList;

		struct RelationNode {
			int m_curSolId = -1;
			RelationNode* m_after = nullptr;
			RelationNode* m_before = nullptr;
			RelationList* m_list = nullptr;

			void init() {
				m_curSolId = -1;
				m_after = nullptr;
				m_before = nullptr;
			}

			void set(int curId) {
				m_curSolId = curId;
			}
			void removeFromLink();
			void insertToLink(RelationNode* head) {
				removeFromLink();
				m_after = head->m_after;
				m_before = head;
				m_after->m_before = this;
				head->m_after = this;
			}

		};


		struct RelationList {
			std::unique_ptr<RelationNode> m_head;
			std::unique_ptr<RelationNode> m_tail;
			int m_num = 0;
			//double m_sumDis2parent = 0;

			void initialize() {
				m_head.reset(new RelationNode);
				m_head.reset(new RelationNode);
				m_head->m_after = m_tail.get();
				m_tail->m_before = m_head.get();
				m_num = 0;
			//	m_sumDis2parent = 0;
			}

			void insertNode(RelationNode* cur);
			void removeNode(RelationNode* cur);

		};


		class NBN_indi : public eax_tsp::TIndi {
		protected:

		public:

			int m_solId = -1;
			std::unique_ptr<RelationNode> m_parentNode;
			RelationList m_sons;

			double m_dis2parent = std::numeric_limits<double>::max();
			double m_fitness = 0;
			int m_stagnation_time = 0;
			bool m_improve = false;
			bool m_stag = false;
			

			double m_maxRadius = 0;

			//unsigned long long m_curVisited = 0;
			unsigned long long m_curClassFlag = 0;


			~NBN_indi() = default;

			void initialize(int N, int solId) {
				define(N);
				m_solId = solId;
				m_parentNode.reset(new RelationNode);
				m_parentNode->init();
				m_dis2parent = std::numeric_limits<double>::max();
				m_fitness = 0;
				m_stagnation_time = 0;
				m_curClassFlag = 0;

				m_sons.initialize();
			}

			void setSolId(int solId) {
				m_solId = solId;
			}

			void copySol(const NBN_indi& other) {
				fLink = other.fLink;
				fEvaluationValue = other.fEvaluationValue;
			}

			void clearParent() {
				m_parentNode->removeFromLink();
				m_parentNode->m_curSolId = -1;
			}
			
			void updateParent(NBN_indi& parent, double minDis) {
				m_parentNode->set(parent.m_solId);
				parent.m_sons.insertNode(m_parentNode.get());
				m_dis2parent = minDis;
				m_stagnation_time = 0;
			}
			void updateParent(NBN_indi& parent) {
				m_parentNode->set(parent.m_solId);
				parent.m_sons.insertNode(m_parentNode.get());
				//m_dis2parent = minDis;
				m_stagnation_time = 0;
			}

			
			int judgeInside(NBN_indi& parent, double dis) {
				if (parent.m_fitness > m_fitness) {
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
				const NBN_indi& a, 
				const NBN_indi& b,
				double& dis
				) {
				dis = -1;
				if (a.m_fitness < b.m_fitness) {
					dis = a.distanceTo(b);
					if (dis == a.m_dis2parent) {
						return -1;
					}
					else if (dis < a.m_dis2parent) {
						return -2;
					}
				}
				else if (a.m_fitness > b.m_fitness) {
					dis = a.distanceTo(b);
					if (dis == b.m_dis2parent) {
						return 1;
					}
					else if (dis < b.m_dis2parent) {
						return 2;
					}
				}
				return 0;
			}

			static int judgeRelationshipWithDis(
				const NBN_indi& a,
				const NBN_indi& b,
				double dis
			) {
				if (a.m_fitness < b.m_fitness) {
					if (dis == a.m_dis2parent) {
						return -1;
					}
					else if (dis < a.m_dis2parent) {
						return -2;
					}
				}
				else if (a.m_fitness > b.m_fitness) {
					if (dis == b.m_dis2parent) {
						return 1;
					}
					else if (dis < b.m_dis2parent) {
						return 2;
					}
				}
				return 0;
			}


			static int judgeRelationshipFit(
				const NBN_indi& better,
				const NBN_indi& worse,
				double dis
			) {
				{
					if (dis == worse.m_dis2parent) {
						return 1;
					}
					else if (dis < worse.m_dis2parent) {
						return 2;
					}
				}
				return 0;
			}

		

			static void updateRelationship(
				NBN_indi& a,
				NBN_indi& b,
				int rf, double curdis,
				ofec::Random* rnd) {
				if (rf == 1) {
					if (rnd->uniform.next() < 0.5) {
						b.updateParent(a, curdis);
					}
				}
				else if (rf == 2) {
					b.updateParent(a, curdis);
				}
				else if (rf == -1) {
					if (rnd->uniform.next() < 0.5) {
						a.updateParent(b, curdis);
					}
				}
				else if (rf == -2) {
					a.updateParent(b, curdis);
				}
			}




		};

	
	public:
		PopNBN_EAX_V2() : Population(), TEnvironment() {}
		PopNBN_EAX_V2(size_t size_pop, Problem* pro) : Population(size_pop, pro), TEnvironment() {};

		virtual ~PopNBN_EAX_V2() = default;
		int evolve(Problem* pro, Algorithm* alg, Random* rnd) override;
		virtual void initialize(Problem* pro, Random* rnd) override;
		//	int evaluate(Problem* pro, Algorithm* alg) override;

		void increaseClassFlag() {
			if (++m_curClassFlag == m_maxFlag) {
				for (auto& it : m_id2inds) {
					it->m_curClassFlag = 0;
				}
				m_curClassFlag = 1;
			}
		}
		//void increaseVisitedFlag() {
		//	if (++m_curVisited == m_maxFlag) {
		//		for (auto& it : m_id2inds) {
		//			it.m_curVisited = 0;
		//		}
		//		m_curVisited = 1;
		//	}
		//}
		

		void updateRelationShipOfV(
			std::vector<int>& newIndi, ofec::Random* rnd) {
			double curdis(0);
			int rf(0);
			for (int idx(0); idx < newIndi.size(); ++idx) {
				auto a = newIndi[idx];
				for (int idy(idx + 1); idy < newIndi.size(); ++idy) {
					auto b = newIndi[idy];
					rf = NBN_indi::judgeRelationship(*m_id2inds[a], *m_id2inds[b], curdis);
					NBN_indi::updateRelationship(*m_id2inds[a], *m_id2inds[b], rf, curdis, rnd);

				}
			}
		}

		
		void updateRelationShipOf2V(
			std::vector<int>& newIndi,
			std::vector<int>& oldIndi,
			ofec::Random* rnd) {
			double curdis(0);
			int rf(0);
			for (int idx(0); idx < newIndi.size(); ++idx) {
				auto a = newIndi[idx];
				for (int idy(0); idy < oldIndi.size(); ++idy) {
					auto b = oldIndi[idy];
					rf = NBN_indi::judgeRelationship(*m_id2inds[a], *m_id2inds[b], curdis);
					NBN_indi::updateRelationship(*m_id2inds[a], *m_id2inds[b], rf, curdis, rnd);

				}
			}
		}
	protected:
		unsigned long long m_maxNumStagation = 20;
		unsigned long long m_curIter = 0;
		GL_calculator m_gl_calculator;
		double m_minDis = 5;
		std::vector<int> m_curPop;
		std::vector<int> m_stagnationIndis;
		std::vector<int> m_history_indis;
		std::vector<int> m_peaks;

		std::vector<std::unique_ptr<NBN_indi>> m_id2inds;
		std::vector<int> m_id2curId;
		
		
		int m_maxPop = 100;
		int m_judgeStagnationIter = 1;
		int m_judgeIter = 5;
		
		int m_curSolId = 0;

		int m_maxStag = 30;
		

		unsigned long long m_curClassFlag = 0;
	//	unsigned long long m_curVisited = 0;
		unsigned long long m_maxFlag = std::numeric_limits<unsigned long long>::max();




		std::vector<eax_tsp::TIndi> m_bestSols; 
		double m_bestObj = 0;

		std::function<void (NBN_indi& sol, Problem* pro)> m_fitnessUpdate;
		


	protected:

		

	public:

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