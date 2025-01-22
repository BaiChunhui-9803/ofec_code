#include "eax_tsp_basin_pop.h"
#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include <thread>


void ofec::PopEaxTspBasin::setNumberThreads(int numThread, Random* rnd) {
	m_threadInfos.resize(numThread);

	if (numThread > m_numThread) {
		for (int idx(m_numThread); idx < numThread; ++idx) {
			initThreadInfo(m_threadInfos[idx],rnd);
		}
	}
	m_numThread = numThread;
}

void ofec::PopEaxTspBasin::initThreadInfo(ThreadInfo& cur, Random* rnd) {
	using namespace ofec::eax_tsp_mt;
	
	auto& N = fEvaluator->Ncity;
	cur.tCross.reset(new TCross(N));
	cur.tCross->eval = fEvaluator;
	cur.tCross->Npop = Npop;

	cur.tKopt.reset(new TKopt(N));
	cur.tKopt->eval = fEvaluator;
	cur.tKopt->tRand = cur.tCross->tRand;
	cur.tKopt->setInvNearList();

	cur.m_random.reset(new Random(rnd->uniform.next()));
}

void ofec::PopEaxTspBasin::initializeSingleThread(const HNSWbasin& model, Environment* env, Random* rnd)
{
	using namespace eax_tsp_mt;
	m_terminatedFlag = false;



	m_pos = (env->problem()->optimizeMode(0) == OptimizeMode::kMaximize) ? 1 : -1;
	m_curIter = 0;
	Npop = 100;
	Nch = 30;

	define(env, rnd->getSharedPtr());
	init();

	setNumberThreads(std::thread::hardware_concurrency(), rnd);
	
	generateRandomSolutions(model, env, rnd);
	selectSolutionsFromBasin(model, env, rnd);
	updateEdgeFrequency();
	setAverageBest(m_curPop);
	
	m_terminatedFlag = eax_tsp_mt::TEnvironment::terminationCondition(m_curPop);
}


void ofec::PopEaxTspBasin::initializeMultiThread(const HNSWbasin& model, Environment* env, Random* rnd)
{
	using namespace eax_tsp_mt;
	m_terminatedFlag = false;

	m_pos = (env->problem()->optimizeMode(0) == OptimizeMode::kMaximize) ? 1 : -1;
	m_curIter = 0;
	Npop = 100;
	Nch = 30;
	define(env, rnd->getSharedPtr());
	init();


	setNumberThreads(std::thread::hardware_concurrency(),rnd);
	generateRandomSolutionMultiThread(model, env, rnd);
	selectSolutionsFromBasin(model, env, rnd);
	updateEdgeFrequency();

	setAverageBest(m_curPop);

	m_terminatedFlag = eax_tsp_mt::TEnvironment::terminationCondition(m_curPop);
	//generateRandomSolutions(model, env, rnd);

}

bool ofec::PopEaxTspBasin::judgeInsideBasin(const HNSWbasin& model, SolutionBase& sol, Environment* env) const{
	return model.calculateBasinId(sol) == m_curBasinId;
}
void ofec::PopEaxTspBasin::generateSolutionsInBasinCross(
	const std::vector<int>& shuffleIds,
	std::vector<std::shared_ptr<SolutionType>>& sonSols,
	std::vector<bool>& activeFlag,
	ThreadInfo& thrdInfo,
	std::mutex& mtx,
	int& maxSolId,
	const HNSWbasin& model, Environment* env) {
	auto tspPro = CAST_TSP(env->problem());
	int curSolId = -1;

	std::vector<std::vector<int>> tmp_fEdgeFreq = fEdgeFreq;
	while (true) {
		{
			std::unique_lock lock(mtx);
			if (maxSolId >= sonSols.size())return;
			curSolId = maxSolId++;
		}

		auto& cursol = sonSols[curSolId];
		{
			//m_id2inds.push_back(NBN_indi());
			cursol.reset(new SolutionType);
			cursol->define(tspPro->numberVariables());
			cursol->copy(*m_curPop[curSolId]);

		}
		 {
			thrdInfo.tCross->setParents(*cursol, *m_curPop[shuffleIds[curSolId]], fFlagC, Nch, thrdInfo.m_random.get());
			thrdInfo.tCross->doIt(*cursol, *m_curPop[shuffleIds[curSolId]], Nch, 1, fFlagC, tmp_fEdgeFreq, thrdInfo.m_random.get());
			//updateFitness(sonSols[idx].get());
			cursol->updateSolBase();
			updateFitness(cursol.get());
		}

		 if(cursol->fitness()>m_curPop[curSolId]->fitness())
		 activeFlag[curSolId] = judgeInsideBasin(model, *cursol, env);

	}
	
}
void ofec::PopEaxTspBasin::generateRandomSolutionOneThread(
	ThreadInfo& thrdInfo,
	std::mutex& mtx,
	int& maxLoop,
	const HNSWbasin& model, Environment* env) {
	auto tspPro = CAST_TSP(env->problem());
	std::shared_ptr<SolutionType> cursol;
	bool isInside = false;
	



	while (true) {
		{
			std::unique_lock lock(mtx);
			if (maxLoop <= 0 || m_curPop.size() >= Npop) return;
			--maxLoop;
		}

		cursol.reset(new SolutionType());
		cursol->define(tspPro->numberVariables());

		thrdInfo.tKopt->makeRandSol(*cursol, thrdInfo.m_random.get());/* Make a random tour */
		thrdInfo.tKopt->doIt(*cursol, thrdInfo.m_random.get()); /* Apply the local search with the 2-opt neighborhood */

		//tKopt->makeRandSol(*cursol); /* Make a random tour */
		//tKopt->doIt(*cursol);        /* Apply the local search with the 2-opt neighborhood */

		updateFitness(cursol.get());
		cursol->updateSolBase();
		isInside = judgeInsideBasin(model, *cursol, env);

		if(isInside){
			std::unique_lock lock(mtx);
			if (m_curPop.size() < Npop) {
				m_history_data.emplace_back(std::move(cursol));
				m_curPop.push_back(m_history_data.back().get());
			}
		}
	}
}
void ofec::PopEaxTspBasin::generateRandomSolutionMultiThread(
	const HNSWbasin& model, Environment* env, Random* rnd) {

	std::mutex mtx;
	int maxLoop = m_maxLoopInGeneration;
	std::vector<std::thread> thrds;
	int num_task = m_numThread;
	//	num_task = 1;
	for (size_t i = 0; i < num_task; ++i) {
		thrds.push_back(std::thread(
			&PopEaxTspBasin::generateRandomSolutionOneThread, 
			this, std::ref(m_threadInfos[i]), std::ref(mtx), std::ref(maxLoop), std::ref(model), env));
	}
	for (auto& thrd : thrds)
		thrd.join();
}
void ofec::PopEaxTspBasin::judgeInsideBasinOneThread(
	std::vector<SolutionType*>& sols, std::vector<bool>& flagInside,
	std::mutex& mtx, int& maxSolId,
	const HNSWbasin& model, Environment* env, Random* rnd
) {
	
	int curSolId = -1;
	while (true) {
		{
			std::unique_lock lock(mtx);
			if (maxSolId >= sols.size())return;
			curSolId = maxSolId++;
		}

		auto& cursol = *sols[curSolId];
		flagInside[curSolId] = judgeInsideBasin(model, cursol, env);
		
	}
}

void ofec::PopEaxTspBasin::judgeInsideBasinMultiThread(
	std::vector<SolutionType*>& sols, std::vector<bool>& flagInside,
	const HNSWbasin& model, Environment* env, Random* rnd
) {

	std::mutex mtx;
	int maxSolId = 0;
	std::vector<std::thread> thrds;
	int num_task = m_numThread;
	//	num_task = 1;
	for (size_t i = 0; i < num_task; ++i) {
		thrds.push_back(std::thread(
			&PopEaxTspBasin::judgeInsideBasinOneThread,
			this, std::ref(sols), std::ref(flagInside), std::ref(mtx), std::ref(maxSolId),  std::ref(model), env, rnd));
	}
	for (auto& thrd : thrds)
		thrd.join();
}

void ofec::PopEaxTspBasin::generateRandomSolutions(const HNSWbasin& model, Environment* env, Random* rnd) {

	auto tspPro = CAST_TSP(env->problem());

	std::shared_ptr<SolutionType> cursol;
	int loop = m_maxLoopInGeneration;
	while (m_curPop.size() < Npop&&--loop) {
		cursol.reset(new SolutionType());
		cursol->define(tspPro->numberVariables());

		tKopt->makeRandSol(*cursol, rnd); /* Make a random tour */
		tKopt->doIt(*cursol, rnd);        /* Apply the local search with the 2-opt neighborhood */
		
		updateFitness(cursol.get());
		cursol->updateSolBase();
		if (judgeInsideBasin(model, *cursol,env)) {
			m_history_data.emplace_back(std::move(cursol));
			m_curPop.push_back(m_history_data.back().get());
		}
	}
}

void ofec::PopEaxTspBasin::updateEdgeFrequency() {
	std::vector<std::vector<std::vector<int>>*> links;
	for (auto& it : m_curPop) {
		links.push_back(&it->fLink);
	}
	calEdgeFreq(links);
}
void ofec::PopEaxTspBasin::selectSolutionsFromBasin(const HNSWbasin& model, Environment* env, Random* rnd) {

	auto tspPro = CAST_TSP(env->problem());
	std::shared_ptr<SolutionType> cursol;

	if (m_curPop.size() < Npop) {
		std::vector<SolutionBase*> candidates;
		for (int idx(0); idx < model.belongBasinId().size(); ++idx) {
			if (model.belongBasinId()[idx] == m_curBasinId) {
				candidates.push_back(model.hnswModel().solutions()[idx].get());
			}
		}
		rnd->uniform.shuffle(candidates.begin(), candidates.end());

		while (m_curPop.size() < Npop) {
			cursol.reset(new SolutionType());
			cursol->define(tspPro->numberVariables());
			cursol->copy(*candidates.back());
			m_history_data.push_back(std::move(cursol));
			m_curPop.push_back(m_history_data.back().get());
			candidates.pop_back();
		}
	}
}
void ofec::PopEaxTspBasin::evolveSingleThread(const HNSWbasin& model, Environment* env, Random* rnd) {
	auto tspPro = CAST_TSP(env->problem());


	// generete solutions
	std::vector<int> shuffleIds(m_curPop.size());
	for (int idx(0); idx < shuffleIds.size(); ++idx) {
		shuffleIds[idx] = idx;
	}
	rnd->uniform.shuffle(shuffleIds.begin(), shuffleIds.end());
	std::vector<std::shared_ptr<SolutionType>> sonSols(m_curPop.size());

	for (int idx(0); idx < sonSols.size(); ++idx) {
		//m_id2inds.push_back(NBN_indi());
		sonSols[idx].reset(new SolutionType);
		sonSols[idx]->define(tspPro->numberVariables());
		sonSols[idx]->copy(*m_curPop[idx]);

	}
	for (int idx(0); idx < sonSols.size(); ++idx) {
		tCross->setParents(*sonSols[idx], *m_curPop[shuffleIds[idx]], fFlagC, Nch,rnd);
		tCross->doIt(*sonSols[idx], *m_curPop[shuffleIds[idx]], Nch, 1, fFlagC, fEdgeFreq,rnd);
		//updateFitness(sonSols[idx].get());
		sonSols[idx]->updateSolBase();
		updateFitness(sonSols[idx].get());
	}

	std::vector<bool> activeFlag(sonSols.size(), false);

	//{
	//	std::vector<int> toSolId(sonSols.size(), -1);
	//	std::vector< SolutionType*> sols;
	//	int fromId = 0;
	//	for (int idx(0); idx < sonSols.size(); ++idx) {
	//		if (sonSols[idx]->fitness() > m_curPop[idx]->fitness()) {
	//			toSolId[idx] = fromId++;
	//			sols.push_back(sonSols[idx].get());
	//		}
	//	}
	//	std::vector<bool> subActiveFlag(fromId, false);
	//	judgeInsideBasinMultiThread(sols, subActiveFlag, model, env, rnd);
	//	
	//	for (int idx(0); idx < sonSols.size(); ++idx) {
	//		if (toSolId[idx] >= 0) {
	//			activeFlag[idx] = subActiveFlag[toSolId[idx]];
	//		}
	//	}
	//}
	
	for (int idx(0); idx < sonSols.size(); ++idx) {
		if (judgeInsideBasin(model, *sonSols[idx], env)&& sonSols[idx]->fitness() > m_curPop[idx]->fitness()) {
			activeFlag[idx] = true;
		}
	}

	for (int idx(0); idx < sonSols.size(); ++idx) {
		if (activeFlag[idx]) {
			m_history_data.push_back(std::move(sonSols[idx]));
			m_curPop[idx] = m_history_data.back().get();
		}
	}

	updateEdgeFrequency();
	setAverageBest(m_curPop);
	m_terminatedFlag = eax_tsp_mt::TEnvironment::terminationCondition(m_curPop);
}


void ofec::PopEaxTspBasin::evolveMultiThread(const HNSWbasin& model, Environment* env, Random* rnd){

	// generete solutions
	std::vector<int> shuffleIds(m_curPop.size());
	for (int idx(0); idx < shuffleIds.size(); ++idx) {
		shuffleIds[idx] = idx;
	}
	rnd->uniform.shuffle(shuffleIds.begin(), shuffleIds.end());
	std::vector<std::shared_ptr<SolutionType>> sonSols(m_curPop.size());
	std::vector<bool> activeFlag(m_curPop.size(), false);

	{
		std::mutex mtx;
		int maxSolId = 0;
		std::vector<std::thread> thrds;
		int num_task = m_numThread;
		//	num_task = 1;
		for (size_t i = 0; i < num_task; ++i) {
			thrds.push_back(std::thread(
				&PopEaxTspBasin::generateSolutionsInBasinCross,
				this, std::cref(shuffleIds), std::ref(sonSols), std::ref(activeFlag), std::ref(m_threadInfos[i]), 
				std::ref(mtx), std::ref(maxSolId), std::cref(model), env));
		}
		for (auto& thrd : thrds)
			thrd.join();
	}
	

	for (int idx(0); idx < sonSols.size(); ++idx) {
		if (activeFlag[idx]) {
			m_history_data.push_back(std::move(sonSols[idx]));
			m_curPop[idx] = m_history_data.back().get();
		}
	}

	updateEdgeFrequency();
	setAverageBest(m_curPop);

	m_terminatedFlag = eax_tsp_mt::TEnvironment::terminationCondition(m_curPop);
}