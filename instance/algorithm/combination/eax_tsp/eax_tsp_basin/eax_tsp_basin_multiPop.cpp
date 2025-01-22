#include "eax_tsp_basin_multipop.h"
#include "../../../../problem/combination/basin_divisioin/init_pop_basin.h"
#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS


void ofec::EaxTspBasinMultiPop::initialize_(Environment* env) {

	Algorithm::initialize_(env);
	auto basinPro = CAST_NBN_Basin(env->problem());
	///m_hnswBasin.initialize(basinPro->hnswModel(), basinPro->belongBasinIds());
	m_hnswBasin.initialize(basinPro->hnswModel(), basinPro->belongBasinIds(), basinPro->basinIds());
	//m_hnswModel.copy(basinPro->hnswModel());
	size_t numPop = m_hnswBasin.basinIds().size();
	m_pops.resize(numPop);
	for (int idx(0); idx < numPop; ++idx) {
		m_pops[idx].setBasinId(idx);
		m_pops[idx].initializeMultiThread(m_hnswBasin, env, m_random.get());
		m_pops[idx].setMaxStagnationIter(m_maxStagnation);
	}

}


void ofec::EaxTspBasinMultiPop::udpateBestSolution() {
	if (m_pops.size() == 0)return;
	for (int idx(0); idx < m_pops.size(); ++idx) {
		auto& curpop = m_pops[idx];
		curpop.updateBestSolution();
		if (curpop.bestSolution() != nullptr) {
			if (m_best == nullptr) {
				m_best = curpop.bestSolution()->getSharedPtr();
			}
			else if (m_best->fitness() < curpop.bestSolution()->fitness()) {
				m_best = curpop.bestSolution()->getSharedPtr();
			}
		}
	}
	//for (auto it = m_pops.begin(); it != m_pops.end(); ++it) {
	//	(*it)->updateBestSolution();
	//	if ((*it)->bestSolution() != nullptr) {
	//		if (m_best == nullptr) {
	//			m_best = (*it)->bestSolution()->getSharedPtr();
	//		}
	//		else if (m_best->fitness() < (*it)->bestSolution()->fitness()) {
	//			m_best = (*it)->bestSolution()->getSharedPtr();
	//		}
	//	}

	//}


}
bool ofec::EaxTspBasinMultiPop::terminating() {
	bool flagNotTerminate = false;
	for (auto& it : m_pops) {
		if (!it->terminationCondition()) {
			flagNotTerminate = true;
			break;
		}
	}
	return !flagNotTerminate;
	//return Algorithm::terminating()|| m_pop->terminationCondition();
	//return Algorithm::terminating() || m_pop->terminationCondition();
}
void ofec::EaxTspBasinMultiPop::run_(Environment* env) {
	while (!terminating()) {

		for (auto it = m_pops.begin(); it != m_pops.end(); ++it) {
			if (!(*it)->terminationCondition()) {
				(*it)->evolveMultiThread(m_hnswBasin, env, m_random.get());
			}
		}
		//m_pop->evolveMultiThread(m_hnswBasin, env, m_random.get());

#ifdef OFEC_DATUM_MULTI_POP_H
		datum::g_multi_pop.resize(m_pops.size());
		for (int idx(0); idx < datum::g_multi_pop.size(); ++idx) {
			auto& curpop = m_pops[idx].curPop();
			auto& datumPop = datum::g_multi_pop[idx];
			datumPop.clear();
			for (auto& it : curpop) {
				datumPop.push_back(it);
	        }
		}
#endif


		datumUpdated(env);
	}

}
