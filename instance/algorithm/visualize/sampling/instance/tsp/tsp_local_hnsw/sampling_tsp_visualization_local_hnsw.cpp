#include "sampling_tsp_visualization_local_hnsw.h"
#include <thread>

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS



void ofec::SamplingTSP_Visualization_localHnsw::initialize_(Environment* env)
{
	SamplingTSP_Visualization::initialize_(env);

	auto pro = CAST_NBN_Basin(env->problem());
	
	m_model.initialize(env, m_random.get(), std::thread::hardware_concurrency());
	if (pro != nullptr) {
		auto& curSolIds = m_solIds[m_numRun];

		std::set<int> filteredSols;
		for (auto& it : curSolIds) {
			for (auto& it2 : it) {
				filteredSols.insert(it2);
			}
		}
		
		auto& totalSols = pro->hnswModel().solutions();
		std::vector<SolutionBase*> sols;
		for (auto& it : filteredSols) {
		//	for (auto& it2 : it) {
				sols.push_back(totalSols[it].get());
		//	}
		}

		//std::cout << "total sol\t" << sols.size() << std::endl;
		std::vector<int> solIds;
		m_model.addDataMutliThread(sols, solIds);
		

		auto& vFitness(m_fitness);
		vFitness.clear();
		//std::vector<Real> vFitness;
		for (auto& it : sols) {
			vFitness.push_back(it->fitness());
		}

		{
			nbn::HnswModelNBN model2;
			model2.copy(m_model);
			//model2.setNumberThreads(1);
			model2.updateFitness(vFitness);
			model2.calculateNBN(m_belong, m_dis2parent);
		}
		m_localNBN = true;
	}


#ifdef OFEC_DATUM_OFFLINE_LOCAL_NBN_H
	ofec::datum::m_algorithm_nbn.m_working = true;
	ofec::datum::m_algorithm_nbn.m_fitness = m_fitness;
	ofec::datum::m_algorithm_nbn.m_belong = m_belong;
	ofec::datum::m_algorithm_nbn.m_dis2parent = m_dis2parent;
	ofec::datum::m_algorithm_nbn.m_sols = m_model.solutions();
	
#endif // OFEC_DATUM_OFFLINE_LOCAL_NBN_H



	
}

void ofec::SamplingTSP_Visualization_localHnsw::run_(Environment* env)
{
	SamplingTSP_Visualization::run_(env);
}

void ofec::SamplingTSP_Visualization_localHnsw::addInputParameters()
{
}
