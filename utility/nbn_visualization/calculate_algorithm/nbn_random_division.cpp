#include "nbn_random_division.h"
#include "../function/custom_function.h"
#include "nbn_grid_division.h"
namespace ofec {


	//void NBN_record::addNode(int id,
	//	const std::shared_ptr<SolBase>& sol,
	//	double fitness, int belong_id, double dis2parent, bool flag_opt 
	//) {

	////	m_list.em

	//	
	//	node cur;
	//	cur.m_id = id;
	//	cur.m_sol = sol;
	//	cur.m_fitness = fitness;
	//	cur.m_belong_id = belong_id;
	//	cur.m_dis2parent = dis2parent;
	//	cur.m_flag_opt = flag_opt;
	//	
	//	m_list.push_back(cur);
	//}

	void NBN_RandomDivision::addRandomSol(const SolBase& originSol, int& id, bool flag_opt,
		int popIter, int popSolId, int algId) {
	//	for(auto& solIter:m_sol)


		// add Node 
		std::shared_ptr<SolBase> sol(m_pro->createSolution(originSol));


		int cur_id(m_list.size());
		double dis2parent = std::numeric_limits<double>::max();

		double minDis = std::numeric_limits<double>::max();
		int belong_id = cur_id;
		double fitness = sol->fitness();
		id = cur_id;
		for (auto& it : m_list) {
			minDis = std::min(minDis, m_pro->norVariableDistance(*it.m_sol, *sol));
		//	if (flag_opt && it.m_flag_opt) continue;
			//if (flag_opt || it.m_fitness < fitness) {
			//	double dis = m_problem->norVariableDistance(*it.m_sol, *sol);
			//	//double dis = it.m_sol->norvar(*sol, m_id_pro);
			//	if (it.m_dis2parent > dis || (it.m_dis2parent == dis && m_rand->uniform.next() < 0.5)) {
			//		it.m_dis2parent = dis;
			//		it.m_belong_id = cur_id;
			//	}
			//}
			//else if (it.m_flag_opt || it.m_fitness > fitness) {
			//	double dis = m_problem->norVariableDistance(*it.m_sol, *sol);

			//	//double dis = it.m_sol->variableDistance(*sol, m_id_pro);
			//	if (dis2parent > dis || (dis2parent == dis && m_rand->uniform.next() < 0.5)) {
			//		dis2parent = dis;
			//		belong_id = it.m_id;
			//	}
			//}


			if ( it.m_fitness < fitness) {
				double dis = m_pro->norVariableDistance(*it.m_sol, *sol);
				//double dis = it.m_sol->norvar(*sol, m_id_pro);
				if (it.m_dis2parent > dis || (it.m_dis2parent == dis && m_random->uniform.next() < 0.5)) {
					it.m_dis2parent = dis;
					it.m_belong_id = cur_id;
				}
			}
			else if (it.m_fitness > fitness) {
				double dis = m_pro->norVariableDistance(*it.m_sol, *sol);

				//double dis = it.m_sol->variableDistance(*sol, m_id_pro);
				if (dis2parent > dis || (dis2parent == dis && m_random->uniform.next() < 0.5)) {
					dis2parent = dis;
					belong_id = it.m_id;
				}
			}
		}


		if (minDis < 1e-6 )return;

		m_list.emplace_back(cur_id, sol, fitness, belong_id,
			dis2parent, flag_opt, popIter, popSolId, algId);
		//addNode(cur_id, sol, fitness, belong_id, dis2parent, flag_opt);

	}

	void NBN_RandomDivision::initialize_(bool flag_grid_sample ) {
		if (flag_grid_sample && m_pro->hasTag(ofec::ProTag::kConOP)) {
			NBN_GridDivision sample_alg;
			sample_alg.setMaxSampleSize(m_maxSample);
			sample_alg.setMultiThread(m_flag_multiThread);
			sample_alg.initialize(m_pro, m_random, m_eval_fun);

			//sample_alg.init(m_problem.get(), m_random.get(), m_eval_fun);
			//sample_alg.initNodesNormal();
			//sample_alg.updateFitness();
			//sample_alg.calculate();
			//sample_alg.udpate_network();

			std::vector<std::shared_ptr<SolBase>> sols;
			std::vector<std::shared_ptr<SolBase>> representatives;
			std::vector<double> fitness;
			std::vector<int> belong;
			std::vector<double> dis2par;
			sample_alg.getNearestBetterNetworkShareMemory(sols, representatives, fitness, belong, dis2par);
			for (int idx(0); idx < sols.size(); ++idx) {
				// addNode
				m_list.emplace_back(idx, sols[idx], fitness[idx], belong[idx], dis2par[idx]);
			}
		}
		else {
			addRandomSols(m_maxSample);
		}
	}



	void NBN_RandomDivision::addRandomSols(int initNum) {
		if (!m_flag_multiThread || m_pro->hasTag(ProTag::kCSIWDN)) {
			SolBase* ptr_cur_sol;
			int id(0);
			for (int idx(0); idx < initNum; ++idx) {
				ptr_cur_sol = m_pro->createSolution();
				m_pro->initSolution(*ptr_cur_sol, m_random.get());
				m_eval_fun(*ptr_cur_sol, m_pro.get());
				std::shared_ptr<SolBase> cur_sol;
				cur_sol.reset(ptr_cur_sol);
				addRandomSol(*cur_sol, id);
			}
		}
		else {
			int id(0);
			std::vector<SolBase*> sols(initNum);
			UTILITY::generateRandomSolutionsMultiThreads(sols, m_random.get(), m_pro.get(), m_eval_fun);
			for (int idx(0); idx < sols.size(); ++idx) {
				std::shared_ptr<SolBase> cur_sol;
				cur_sol.reset(sols[idx]);
				addRandomSol(*cur_sol, id);
			}
			sols.clear();
		}

	}
	void NBN_RandomDivision::getSharedData(
		std::vector<std::shared_ptr<SolBase>>& sols,
		std::vector<double>& fitness, std::vector<int>& belongs)const
	{
		sols.resize(m_list.size());
		belongs.resize(m_list.size());
		fitness.resize(m_list.size());
		for (auto& it : m_list) {
			sols[it.m_id] = it.m_sol;
			belongs[it.m_id] = it.m_belong_id;
			fitness[it.m_id] = it.m_fitness;
		}

	}



	void NBN_RandomDivision::getSharedData(
		std::vector<std::shared_ptr<SolBase>>& sols,
		std::vector<double>& fitness,
		std::vector<int>& belongs,
		std::vector<int> popIters,
		std::vector<int> popSolIds,
		std::vector<int> algIds
	)const {
		sols.resize(m_list.size());
		belongs.resize(m_list.size());
		fitness.resize(m_list.size());
		popIters.resize(m_list.size());
		popSolIds.resize(m_list.size());
		algIds.resize(m_list.size());

		for (auto& it : m_list) {
			sols[it.m_id] = it.m_sol;
			belongs[it.m_id] = it.m_belong_id;
			fitness[it.m_id] = it.m_fitness;
			popIters[it.m_id] = it.m_popIter;
			popSolIds[it.m_id] = it.m_popSolId;
			algIds[it.m_id] = it.m_algId;
		}
	}


	void NBN_RandomDivision::getSharedNBN(
		std::vector<std::shared_ptr<SolBase>>& sols,
		std::vector<double>& fitness,
		std::vector<int>& belongs,
		std::vector<double>& dis2parent,
		std::vector<bool>& flagOpt
	)const {
		sols.resize(m_list.size());
		belongs.resize(m_list.size());
		fitness.resize(m_list.size());
		dis2parent.resize(m_list.size());
		flagOpt.resize(m_list.size());
		for (auto& it : m_list) {
			sols[it.m_id] = it.m_sol;
			belongs[it.m_id] = it.m_belong_id;
			fitness[it.m_id] = it.m_fitness;
			dis2parent[it.m_id] = it.m_dis2parent;
			flagOpt[it.m_id] = it.m_flag_opt;
 		}
	}


	void NBN_RandomDivision::getSharedNBN(
		std::vector<std::shared_ptr<SolBase>>& sols,
		std::vector<double>& fitness,
		std::vector<int>& belongs,
		std::vector<double>& dis2parent,
		std::vector<int> popIters,
		std::vector<int> popSolIds,
		std::vector<int> algIds
	)const {


		sols.resize(m_list.size());
		belongs.resize(m_list.size());
		fitness.resize(m_list.size());
		dis2parent.resize(m_list.size());

		popIters.resize(m_list.size());
		popSolIds.resize(m_list.size());
		algIds.resize(m_list.size());

		for (auto& it : m_list) {
			sols[it.m_id] = it.m_sol;
			belongs[it.m_id] = it.m_belong_id;
			fitness[it.m_id] = it.m_fitness;
			dis2parent[it.m_id] = it.m_dis2parent;

			popIters[it.m_id] = it.m_popIter;
			popSolIds[it.m_id] = it.m_popSolId;
			algIds[it.m_id] = it.m_algId;
		}
	}
}