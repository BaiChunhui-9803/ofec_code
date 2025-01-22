#include "eax_tsp_basin.h"
#include "../../../../problem/combination/basin_divisioin/init_pop_basin.h"
#include <thread>
#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS


void ofec::EaxTspBasin::addInputParameters(){}

void ofec::EaxTspBasin::initialize_(Environment* env){
	
	Algorithm::initialize_(env);
	auto basinPro = CAST_NBN_Basin(env->problem());
	m_hnswBasin.initialize(basinPro->hnswModel(), basinPro->belongBasinIds(), basinPro->basinIds());
	//m_hnswModel.copy(basinPro->hnswModel());


	m_pop.reset(new PopEaxTspBasin);
	m_pop->setBasinId(basinPro->currentBasinId());

	m_pop->initializeMultiThread(m_hnswBasin, env, m_random.get());

	m_pop->setNumberThreads(std::thread::hardware_concurrency(), m_random.get());

}

void ofec::EaxTspBasin::run_(Environment* env){
	while (!terminating()) {
		m_pop->evolveMultiThread(m_hnswBasin, env, m_random.get());

#ifdef OFEC_DATUM_MULTI_POP_H
		datum::g_multi_pop.resize(1);
		auto& curpop = m_pop->curPop();
		datum::g_multi_pop.front().clear();
		for (auto& it : curpop) {
			datum::g_multi_pop.front().push_back(it);
		}
#endif

		
		datumUpdated(env);
	}
	
}
