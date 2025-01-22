#include "eax_tsp_alg.h"
#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../../../../../core/environment/environment.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS


void ofec::EAX_TSP::initialize_(Environment *env) {
	using namespace eax_tsp;
	Algorithm::initialize_(env);
	m_env.reset(new TEnvironment());
	m_env->setNPop(100);
	m_env->setNch(30);
	// optimum = -1;
	m_env->setTMax(0);
	// InitURandom(runId);
	m_env->setTeminate(false);

	m_env->define(env, m_random);
	m_env->setMaxGeneration(1e5);
	m_env->setMaxStagnationIter(300);
	//	curInfo.runId = id_run;
	m_env->initialize();

}

void ofec::EAX_TSP::run_(Environment* env) {
	while (!terminating()) {
#ifdef OFEC_DATUM_TSP_OFFLINE_NBN_H
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

#endif

		datumUpdated(env);



		m_env->evolve();
	}
}
