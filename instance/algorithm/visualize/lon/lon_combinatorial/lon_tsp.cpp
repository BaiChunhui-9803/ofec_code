#include "lon_tsp.h"
#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"

void ofec::LocalOptimaNetworkTSP::samplingLKH(ThreadInfo& cur, GlobalInfo& globalInfo, Problem* pro)
{
	using namespace LKH;
	using namespace ofec;

	auto tsp_pro = CAST_TSP(pro);
	auto& filedir = tsp_pro->filePath();
	auto& filename = tsp_pro->fileName();


	LKH::LKHAlg alg;

	alg.readProblem(filedir, filename);
	alg.set2optLocalSearchParamenters();
	alg.assignedMemory();

	TravellingSalesman::HashSolutionMap solMap;
	solMap.initialize(cur.m_rnd.get(), tsp_pro->numberVariables() + 10);

	int KickType;   /* Specifies K for a K-swap-kick */
	KickType = 4;

	int fromSolId(-1), toSolId(-1);
	GainType fromSolCost(0), toSolCost(0);
	//GainType curSolCost = 0;
	unsigned curRunId = 0;
	std::vector<int> sol;
	//auto& alg = cur.m_alg;
	//RecordInfo curInfo;

	while (true) {
		curRunId = 0;
		{
			globalInfo.m_mtx.lock();
			if (globalInfo.m_maxRun) curRunId = globalInfo.m_maxRun--;
			globalInfo.m_mtx.unlock();
		}
		if (curRunId == 0)break;

		lon::TraceEdge curEdge;
		curEdge.run = curRunId;
		if (curEdge.run % 2) {
			alg.setInitialTourAlgorithm(QUICK_BORUVKA);
			fromSolCost = alg.GreedyTour();
		}
		else {
			fromSolCost = alg.generateSolRandom();
		}
		fromSolCost = alg.LinKernighanLocalSearch();
		alg.RecordBetterTour();
		alg.getBetterSolution(sol);

		fromSolId = solMap.getSolId(sol);

		if (cur.m_sols.size() <= fromSolId) {
			cur.m_sols.push_back(sol);
			cur.m_nodeFit.push_back(fromSolCost);
		}


		for (int iter = 0; iter < chainedLKLoop; ++iter) {
			curEdge.iter = iter;
			alg.setCurrentSolutionFromBest();
			alg.KSwapKick(KickType);
			toSolCost = alg.LinKernighanLocalSearch();
			if (fromSolCost >= toSolCost) {
				alg.RecordBetterTour();
				alg.getBetterSolution(sol);
				toSolId = solMap.getSolId(sol);
				if (cur.m_sols.size() <= toSolId) {
					cur.m_sols.push_back(sol);
					cur.m_nodeFit.push_back(toSolCost);
				}

				curEdge.id_from = fromSolId;
				curEdge.id_to = toSolId;
				curEdge.num_kick = 1;

				fromSolId = toSolId;
				fromSolCost = toSolCost;
			}
			else {
				curEdge.id_from = fromSolId;
				curEdge.id_to = fromSolId;
				curEdge.num_kick = 1;
			}
			cur.m_historyEdge.push_back(curEdge);
		}
	}

	alg.freeAll();
	

}

void ofec::LocalOptimaNetworkTSP::insertDatas(std::vector<ThreadInfo>& totalInfo, Problem* pro, Random* rnd) {
	//using namespace lon;
	using namespace ofec;

	auto tsp_pro = CAST_TSP(pro);
	
	TravellingSalesman::HashSolutionMap solMap;
	solMap.initialize(rnd, tsp_pro->numberVariables() + 10);
	
	clear();
	
	for (auto& curInfo : totalInfo) {
		std::vector<int> curIdToNewId(curInfo.m_sols.size());
		for (int idx(0); idx < curIdToNewId.size(); ++idx) {
			curIdToNewId[idx] = solMap.getSolId(curInfo.m_sols[idx]);
			auto& curId = curIdToNewId[idx];
			if (curId >= m_sols.size()) {
				m_sols.push_back(curInfo.m_sols[idx]);
				m_nodeFit.push_back(curInfo.m_nodeFit[idx]);
			}
		}

		lon::TraceEdge curEdge;
		for (auto& oldEdge : curInfo.m_historyEdge) {
			curEdge = oldEdge;
			curEdge.id_from = curIdToNewId[curEdge.id_from];
			curEdge.id_to = curIdToNewId[curEdge.id_to];
			m_historyEdge.push_back(curEdge);
		}
	}
	

}

void ofec::LocalOptimaNetworkTSP::sampleSingleThread(Problem* pro, Random* rnd)
{
	GlobalInfo globalInfo;
	globalInfo.m_maxRun = Runs;
	ThreadInfo curInfo;
	curInfo.initialize(pro,rnd);

	samplingLKH(curInfo, globalInfo, pro);
	
	m_historyEdge = std::move(curInfo.m_historyEdge);
	m_sols = std::move(curInfo.m_sols);
	m_nodeFit = std::move(curInfo.m_nodeFit);
	
	
	curInfo.clear();

}

void ofec::LocalOptimaNetworkTSP::sampleMultiThread(Problem* pro, Random* rnd)
{
	GlobalInfo globalInfo;
	globalInfo.m_maxRun = Runs;

	int numThread = std::thread::hardware_concurrency();
	std::vector<ThreadInfo> threadInfos(numThread);
	for (auto& it : threadInfos) {
		it.initialize(pro, rnd);
	}


	std::vector<std::thread> thrds;
	int num_task = numThread;
	//	num_task = 1;
	for (size_t i = 0; i < num_task; ++i) {
		thrds.push_back(std::thread(
			&LocalOptimaNetworkTSP::samplingLKH, this, std::ref(threadInfos[i]), std::ref(globalInfo), pro ));
	}
	for (auto& thrd : thrds)
		thrd.join();
	


	insertDatas(threadInfos, pro, rnd);
	

	for (auto& it : threadInfos) {
		it.clear();
	}
}
void ofec::LocalOptimaNetworkTSP::trasferToLon(lon::LonInfo& lon ) {
	std::vector<lon::TraceEdge>  edgeHistorylon(m_historyEdge.size());
	for (int idx(0); idx < m_historyEdge.size(); ++idx) {
		edgeHistorylon[idx].count = m_historyEdge[idx].count;
		edgeHistorylon[idx].id_from = m_historyEdge[idx].id_from;
		edgeHistorylon[idx].id_to = m_historyEdge[idx].id_to;
		edgeHistorylon[idx].iter = m_historyEdge[idx].iter;
		edgeHistorylon[idx].num_kick = m_historyEdge[idx].num_kick;
		edgeHistorylon[idx].run = m_historyEdge[idx].run;
	}

//	edgeHistory.clear();
	std::vector<lon::TraceEdge> edges;
	lon::transfer(edgeHistorylon, edges);

	//filterNetowrkByTrace(trace,edges,nodefit,lon);
	std::set<int> funnelSetIds;
	std::set<int> sinkNodes;
	lon::getFunnelIdx(funnelSetIds, edgeHistorylon);
	////getSinkNode(funnelSetIds, edges, nodefit, sinkNodes);
	lon::getSinkNode3(funnelSetIds, edges, m_nodeFit, sinkNodes);
	////bfsLON(sinkNodes, edges,nodefit,lon);

	//filterNetworkRemoveLoop(edges, nodefit, lon);
	lon::filterNetworkRemoveLoopWishSink(sinkNodes, edges, m_nodeFit, lon);
	//bfsNetwork(sinkNodes, edges, nodefit, lon);
	//lon::LonInfo filterLon;
	//lon::filterLON(lon, filterLon, filterRatio);
}

void ofec::LocalOptimaNetworkTSP::ThreadInfo::initialize(Problem* pro,Random* rnd)
{	
	m_rnd.reset(new Random(rnd->uniform.next()));
}

void ofec::LocalOptimaNetworkTSP::ThreadInfo::clear()
{
	//m_alg.freeAll();
}
