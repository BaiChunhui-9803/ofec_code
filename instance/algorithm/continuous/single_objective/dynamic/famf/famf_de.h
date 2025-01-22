/********* Begin Register Information **********
{
	"name": "FAMFDE",
	"identifier": "FAMFDE",
	"problem tags": [ "ConOP", "SOP", "GOP", "MMOP", "DOP" ]
}
*********** End Register Information **********/

#ifndef FAMFDE_H
#define FAMFDE_H

#include "famf.h"
#include "../../../../template/classic/de/Solution.h"
#include "../../../../template/classic/de/population.h"

namespace ofec
{
		using FAMFDE = FAMF<FAMFPop<PopDE<FAMFSolution<IndDE>>>>;
		//using FAMF_PSO = FAMF<FAMFPop<FAMFPop_PSO>>;

	//	class FAMF_PSO : public FAMF<FAMFPop<FAMFPop_PSO>> {

	//	public:
	//		FAMF_PSO(param_map& v) : FAMF<FAMFPop<FAMFPop_PSO>>(v) {
	//			m_C1 = v.at("accelerator1"), m_C2 = v.at("accelerator2");						//accelerators
	//			m_weight = v.at("weight");								// inertia weight	
	//		}
	//		void initialize() override {
	//			for (int i(0); i < m_sub_pop.size(); i++) m_sub_pop[i].set_parameters(m_weight, m_C1, m_C2);
	//			FAMF::initialize();
	//		}
	//		virtual int create_new_swarms(int Solution_num) {
	//			init_Solutions(Solution_num);
	//			//for (int i(0); i < m_sub_pop.size(); i++) {
	//			//	for (int j(0); j < m_sub_pop[i].size(); j++) {
	//			//		m_new_indis.push_back(m_sub_pop[i][j]);
	//			//	}
	//			//}
	//			//for (int i(0); i < m_sub_pop.size(); i++) {
	//			//	m_sub_pop.pop_back();
	//			//	i--;
	//			//}
	//			HSLH<Solution_type> cluster(m_new_indis);
	//			vector<Solution_type> left_indis;

	//			cluster.clustering(min_subpopSize);
	//			for (int i(0); i < cluster.size(); ++i) {
	//				if (cluster[i].size() > m_min_num_inds) {
	//					auto a = new FAMFPop<FAMFPop_PSO>(cluster[i].size(),m_C1,m_C2,m_weight);
	//					m_sub_pop + (a);
	//					auto iter = cluster[i].begin();
	//					for (int j(0); j < cluster[i].size(); ++j, ++iter) {
	//						m_sub_pop[m_sub_pop.size() - 1][j] = *(iter->second);
	//					}
	//					m_sub_pop[m_sub_pop.size() - 1].estimateInitialRadius();
	//				}
	//				else {
	//					auto iter = cluster[i].begin();
	//					for (int j(0); j < cluster[i].size(); ++j, ++iter) {
	//						left_indis.push_back(*(iter->second));
	//					}
	//				}
	//				//m_sub_pop.add(newPop);
	//			}
	//			swap(m_new_indis, left_indis);
	//			return cluster.size();
	//		}

	//		void measureMultiPop(bool) override{
	//			for (int i(0); i < m_sub_pop.size(); ++i) {
	//				for (int j(0); j < m_sub_pop[i].size(); j++) {
	//					m_sub_pop[i][j].evaluate(true, caller::Algorithm);
	//					m_sub_pop[i][j].pbest() = m_sub_pop[i][j];
	//				}
	//				m_sub_pop[i].update_best();
	//			}
	//		}
	//#ifdef OFEC_DEMO
	//		void update_buffer()
	//		{
	//			if (!this->m_initialized)
	//				return;
	//			vector<vector<solution<>*>> pops(this->m_sub_pop.size());
	//			for (size_t k = 0; k < pops.size(); ++k) {
	//				for (size_t i = 0; i < m_sub_pop[k].size(); ++i)
	//					pops[k].emplace_back(&m_sub_pop[k][i].pbest());
	//			}
	//			dynamic_cast<ofec_demo::buffer_cont*>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops);
	//		}
	//#endif

	//	protected:
	//		Real m_C1, m_C2, m_weight;
	//	};
	//}
}
#endif