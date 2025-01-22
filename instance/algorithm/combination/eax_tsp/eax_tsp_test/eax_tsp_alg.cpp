#include "eax_tsp_alg.h"
#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../../../../../core/environment/environment.h"

#ifdef OFEC_PLAYBACK
#include <playback/global.h>
#include "../../../../../buffer_datum/algorithm/multi_pop.h"
#endif



void ofec::EAX_TSP2::initialize_(Environment *env) {
	using namespace eax_tsp2;
	Algorithm::initialize_(env);
	m_env.reset(new TEnvironment());
	m_env->setNPop(100);
	m_env->setNch(30);
	// optimum = -1;
	m_env->setTMax(0);
	// InitURandom(runId);
	m_env->setTeminate(false);

	m_env->define(env, m_random.get());
	m_env->setMaxGeneration(1e5);
	m_env->setMaxStagnationIter(300);
	//	curInfo.runId = id_run;
	m_env->initialize(m_random.get());

#ifdef OFEC_PLAYBACK
	buffer_datum::g_multi_pop.working = true;
#endif
}

void ofec::EAX_TSP2::run_(Environment* env) {
	while (!terminating()) {
#ifdef OFEC_PLAYBACK
		buffer_datum::g_multi_pop.pops.clear();
		std::vector<std::shared_ptr<SolutionBase>> sols(m_env->getNpop());
		std::vector<int> curSol;
		for (int idIndi(0); idIndi < m_env->getNpop(); ++idIndi) {
			m_env->getCurPop()[idIndi].transferSol(curSol);
			auto& it = sols[idIndi];
			it.reset(env->problem()->createSolution());
			auto& tspSol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*it);
			tspSol.variable().vect() = curSol;
			tspSol.evaluate(env);
		}
		buffer_datum::g_multi_pop.pops.resize(1);
		for (auto& it : sols) {
			buffer_datum::g_multi_pop.pops.front().push_back(it.get());
		}
		ofec_playback::g_buffer_manager->appendAlgBuffer(env);
#endif
		m_env->evolve(m_random.get());
	}
}
