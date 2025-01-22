
#include <filesystem>
#include <string>
#include <vector>
#include <iostream>
#include <limits>

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "../utility/hnsw/hnsw_solution_multithread/hnsw_sol_model.h"
#include "../core/global.h"
#include "../core/problem/optima.h"
#include "../core/definition.h"
#include "../core/global.h"
#include "../utility/nbn_visualization/tree_graph_simple.h"
#include "../utility/nbn_visualization/nbn_grid_tree_division.h"
#include "../utility/nbn_visualization/nbn_grid_tree_division_allSize.h"
#include  "../utility/function/custom_function.h"
#include "../utility/function/general_function.h"
#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator.h"
//#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include "interface.h"
//#include "../utility/function/custom_function.h"
#include "../utility/nbn_visualization/tree_graph_simple.h"
#include "../core/algorithm/population.h"
//#include "../utility/nondominated_sorting/filter_sort.h"
#include "../instance/algorithm/visualize/sampling/instance/sampling_eax_tsp.h"
#include "../instance/algorithm/visualize/sampling/sampling_multiThread.h"
#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include <queue>
#include "../core/definition.h"
#include <deque>
#include <chrono>   
//#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator_multiThead.h"

#include "../utility/hnsw/hnsw_nbn/neighbor_info.h"
#include "../utility/hnsw/hnsw_nbn/select_policy.h"
#include "../utility/hnsw/hnsw_nbn/node.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_model.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_nbn.h"
#include "../core/environment/environment.h"
#include "../core/parameter/variants_to_stream.h"

//#include "../instance/algorithm/combination/lkh/lkh_algorithm.h"
#include<filesystem>
#include<queue>
#include <set>

#include "../core/parameter/variants_to_stream.h"
#include "../instance/algorithm/combination/eax_tsp/eax_tsp_test/eax_tsp_alg.h"
#include "../instance/problem/combination/travelling_salesman/tsp_offline_data/travelling_salesman_offline_data.h"



#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>

#include "../instance/algorithm/combination/LKH_origin/INCLUDE/LKH.h"
#include "../instance/algorithm/visualize/lon/lon_combinatorial/lon_tsp.h"
#include "../utility/function/custom_function.h"


using namespace std;





struct RecordInfo {
	int runId = 0;
	long long beforeFit = 0;
	int beforeId = 0;

	long long afterFit = 0;
	int afterId = 0;
};


struct node {
	int nodeId = 0;
	unsigned hashval;
	std::vector<int> sol;
	LKH::GainType solCost;
};



void runTaskTest() {
	using namespace ofec;
	using namespace LKH;

	std::string tspFileDir = g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string tspName = "2202.tsp";

	LKH::LKHAlg alg;

	alg.readProblem(tspFileDir, tspName);
	alg.set2optLocalSearchParamenters();
	alg.assignedMemory();

	int Runs = 1000;
	int chainedLKLoop = 1e4;
	std::vector<node> m_nodes;
	node curNode;
	std::list<RecordInfo> lon_trait;
	RecordInfo curInfo;
	int maxId = 0;
	LKH::GainType curSolCost = 0;
	unsigned Seed = 0;
	int KickType;   /* Specifies K for a K-swap-kick */
	KickType = 4;

	std::vector<int> sol;

	//TravellingSalesman::HashSolutionTree m_hashTree;


	for (int Run = 1; Run <= Runs; Run++) {

		curInfo.runId = Run;
		if (Run % 2) {
			alg.setInitialTourAlgorithm(QUICK_BORUVKA);
			curSolCost = alg.GreedyTour();
		}
		else {
			curSolCost = alg.generateSolRandom();
		}
		curSolCost = alg.LinKernighanLocalSearch();
		alg.RecordBetterTour();
		alg.getBetterSolution(sol);
		unsigned curhash = alg.calHash(sol);


		int curId(0);
		maxId = std::max(curId, maxId);


		curInfo.beforeFit = curSolCost;
		curInfo.beforeId = curId;

		bool flag_exist(false);
		for (int j = 0; j < chainedLKLoop; ++j) {
			alg.setCurrentSolutionFromBest();

			alg.KSwapKick(KickType);
			curSolCost = alg.LinKernighanLocalSearch();
			if (curInfo.beforeFit >= curSolCost) {
				alg.RecordBetterTour();
				alg.getBetterSolution(sol);
				curhash = alg.calHash(sol);
				//flag_exist = lon_set.insertValue(sol, curhash, curId);
				maxId = std::max(maxId, curId);
				curInfo.afterFit = curSolCost;
				curInfo.afterId = curId;
				lon_trait.push_back(curInfo);


				curInfo.beforeFit = curInfo.afterFit;
				curInfo.beforeId = curInfo.afterId;
			}
		}

		alg.SRandom(++Seed);
	}

	alg.freeAll();


}


void runTask(const std::string& saveDir, const std::string& tspfilename) {

	using namespace ofec;

	ofec::ParameterMap params;
	//params["problem name"] = std::string("TSP");
	params["dataFile1"] = tspfilename;
	//params["algorithm name"] = std::string("EAX-TSP-sampling");
	std::string proname = "TSP";


	std::shared_ptr<Environment> env(Environment::create());
	env->recordInputParameters();
	env->initialize();

	env->setProblem(ofec::Factory<ofec::Problem>::produce(proname));
	env->problem()->inputParameters().input(params);
	env->problem()->recordInputParameters();
	env->initializeProblem(0.5);

	auto pro = env->problem();

	auto rnd = std::make_shared<Random>(0.5);

	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));

	LocalOptimaNetworkTSP lon_sampling;
	lon_sampling.sampleMultiThread(pro, rnd.get());

	lon::LonInfo lon;
	lon_sampling.trasferToLon(lon);

	lon::LonInfo filterLon;
	lon::filterLON(lon, filterLon, 0.001);






	std::ofstream out(saveDir + tspfilename + "_lon_001.txt");

	{
		std::stringstream buf;
		UTILITY::vectorToStream(filterLon.m_node_funnel_fit, buf);
		UTILITY::vectorToStream(filterLon.m_node_size, buf);
		UTILITY::vectorToStream(filterLon.m_node_fit, buf);
		UTILITY::vvectorToStream(filterLon.m_graph, buf);

		out << buf.rdbuf();
	}
	out.close();




}





void outputLon(std::stringstream& buf, lon::LonInfo& filterLon) {
	UTILITY::vectorToStream(filterLon.m_node_funnel_idxs, buf);
	UTILITY::vectorToStream(filterLon.m_node_funnel_fit, buf);
	UTILITY::vectorToStream(filterLon.m_node_size, buf);
	UTILITY::vectorToStream(filterLon.m_node_fit, buf);
	UTILITY::vvectorToStream(filterLon.m_graph, buf);
}


void inputLon(std::stringstream& buf, lon::LonInfo& lon) {
	UTILITY::streamToVector(lon.m_node_funnel_idxs, buf);
	UTILITY::streamToVector(lon.m_node_funnel_fit, buf);
	UTILITY::streamToVector(lon.m_node_size, buf);
	UTILITY::streamToVector(lon.m_node_fit, buf);
	UTILITY::streamToVvector(lon.m_graph, buf);
	lon.m_nodeId.resize(lon.m_node_fit.size());
	for (int idx(0); idx < lon.m_nodeId.size(); ++idx) {
		lon.m_nodeId[idx] = idx;
	}
}

void outputLon(const std::string& filename, lon::LonInfo& filterLon) {

	std::ofstream out(filename);
	{
		std::stringstream buf;
		outputLon(buf, filterLon);
		out << buf.rdbuf();
	}
	out.close();
}

void inputLon(const std::string& filename, lon::LonInfo& filterLon) {

	std::ifstream in(filename);
	{
		std::stringstream buf;
		buf << in.rdbuf();
		inputLon(buf, filterLon);
	}
	in.close();
}

void runTask(const std::string& saveDir,
	const std::string& dataDir,
	const std::string& tspfilename) {

	using namespace ofec;

	ofec::ParameterMap params;
	//params["problem name"] = std::string("TSP");
	params["dataFile1"] = tspfilename;
	params["dataDirectory"] = dataDir;
	//params["algorithm name"] = std::string("EAX-TSP-sampling");
	std::string proname = "TSP";


	std::shared_ptr<Environment> env(Environment::create());
	env->recordInputParameters();
	env->initialize();

	env->setProblem(ofec::Factory<ofec::Problem>::produce(proname));
	env->problem()->inputParameters().input(params);
	env->problem()->recordInputParameters();
	env->initializeProblem(0.5);

	auto pro = env->problem();

	auto rnd = std::make_shared<Random>(0.5);

	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));

	LocalOptimaNetworkTSP lon_sampling;
	lon_sampling.sampleMultiThread(pro, rnd.get());

	lon::LonInfo lon;
	lon_sampling.trasferToLon(lon);

	lon::LonInfo filterLon;
	lon::filterLON(lon, filterLon, 0.001);


	outputLon(saveDir + tspname + "_lon_001.txt", filterLon);

	outputLon(saveDir + tspname + "_lon.txt", lon);

	ParameterVariantStream variantStr;
	variantStr << lon_sampling.sols().size();
	for (auto& it : lon_sampling.sols()) {
		variantStr << it;
	}
	//variantStr << lon_sampling.sols();
	variantStr << lon_sampling.nodeFit();

	{
		std::stringstream buf;
		ofec::variants_stream::parameters2StreamMutithread(buf, variantStr);
		std::ofstream out(saveDir + tspname + "_sol.txt");
		out << buf.rdbuf();
		out.close();
	}

}


void inputLonInfo(const std::string& saveDir, const std::string& tspfilename) {

	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));
	lon::LonInfo lon;
	inputLon(saveDir + tspname + "_lon.txt", lon);

	lon::LonInfo filterLon;
	lon::filterLON(lon, filterLon, 0.005);


	outputLon(saveDir + tspname + "_lon_005_test.txt", filterLon);
}



namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run() {

		using namespace ofec;
		using namespace std;

		registerInstance();
		ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
		ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";

		// runTask();


		ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
		// ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";


		std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
		saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/lonData/";

		std::filesystem::create_directories(saveDir);


		std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
		std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";


		runTask(saveDir, dir2, "u574.tsp");

		runTask(saveDir, dir1, "1281.tsp");
		runTask(saveDir, dir1, "6702.tsp");
		runTask(saveDir, dir1, "2202.tsp");
		runTask(saveDir, dir1, "6717.tsp");
		//inputLonInfo(saveDir, "GR666.tsp");

	}


}