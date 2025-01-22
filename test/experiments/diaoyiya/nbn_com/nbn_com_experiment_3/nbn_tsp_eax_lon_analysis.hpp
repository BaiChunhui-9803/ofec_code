
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
#include "../utility/nbn_visualization/visualization/tree_graph_simple.h"
#include  "../utility/function/custom_function.h"
#include "../utility/function/general_function.h"
#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator.h"
//#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include "interface.h"
//#include "../utility/function/custom_function.h"
#include "../core/algorithm/population.h"
//#include "../utility/nondominated_sorting/filter_sort.h"
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

#include "../utility/nbn_visualization/nbn_fla/nbn_fla.h"
#include "../instance/algorithm/combination/eax_tsp/eax_tsp_basin/eax_tsp_basin_multipop.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
#include "../instance/algorithm/combination/LKH_origin/INCLUDE/LKH.h"


#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>

#include "../instance/algorithm/visualize/lon/lon_combinatorial/lon_tsp.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_info.h"


#include<filesystem>
#include<queue>
#include <set>

#include "../core/parameter/variants_to_stream.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla.h"

#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
#include "../instance/algorithm/combination/LKH_origin/INCLUDE/LKH.h"


#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include "../utility/nbn_visualization/instance/travelling_salemans_problem/nbn_modify_tsp.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_info.h"
#include "../utility/nbn_visualization/nbn_fla/tsp_related/nbn_fla_tsp.h"
#include "../utility/nbn_visualization/nbn_calculator/nbn_improve_accuracy_multithread.h"

#include "../utility/matlab/matlab_utility.h"
#include "../utility/nbn_visualization/calculate_algorithm/nbn_iteration_multithread.h"



// lon related network



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




void genenrateTSPenv(
	const std::string& readDir,
	const std::string& tspfilename,
	std::shared_ptr<ofec::Environment>& env) {
	using namespace ofec;
	ofec::ParameterMap params;
	//params["problem name"] = std::string("TSP");
	params["dataFile1"] = tspfilename;
	params["dataDirectory"] = readDir;

	//params["algorithm name"] = std::string("EAX-TSP-sampling");
	std::string proname = "TSP";

	env.reset(Environment::create());
	env->recordInputParameters();
	env->initialize();
	env->setProblem(ofec::Factory<ofec::Problem>::produce(proname));
	env->problem()->inputParameters().input(params);
	env->problem()->recordInputParameters();
	env->initializeProblem(0.5);
}

#include "../utility/nbn_visualization/nbn_calculator/nbn_improve_accuracy_multithread.h"





void udpateNBNinfoIteration(
	ofec::nbn::NBNinfo& m_nbnInfo,
	int solSize,
	ofec::Environment* env, ofec::Random* rnd) {


	std::vector<bool> updateIds;
	int from(m_nbnInfo.m_belong.size() - solSize);
	int to = m_nbnInfo.m_belong.size();
	std::vector<int> sortedIds(m_nbnInfo.m_belong.size());
	for (int idx(0); idx < sortedIds.size(); ++idx) {
		sortedIds[idx] = idx;
	}
	std::sort(sortedIds.begin(), sortedIds.end(), [&](int a, int b) {
		return m_nbnInfo.m_vFitness[a] > m_nbnInfo.m_vFitness[b];
		});

	ofec::nbn::udpateNBNinfoIterationMultiThread(
		m_nbnInfo, updateIds, from, to, sortedIds, env, rnd);

	double improveAccuracy(0);

	for (int idx = from; idx < to; ++idx) {
		if (updateIds[idx]) ++improveAccuracy;
	}
	std::cout << "improveAccuracy\t" << improveAccuracy << std::endl;


}


void runEax(const std::string& readDir,
	const std::string& tspfilename,
	std::vector<std::vector<std::vector<std::shared_ptr<ofec::SolutionBase>>>>& sols, int numRun = 30) {
	using namespace ofec;
	using namespace std;
	using namespace chrono;


	std::string proname = "TSP";
	std::string algName = "EAX-TSP-sampling";

	ofec::ParameterMap params;
	//params["problem name"] = std::string("TSP");
	params["dataFile1"] = tspfilename;
	params["dataDirectory"] = readDir;

	//params["algorithm name"] = std::string("EAX-TSP-sampling");


	//solbases.clear();
	std::vector<std::shared_ptr<ofec::SolutionBase>> tmp;
	std::vector<std::shared_ptr<ofec::SolutionInfo>> solInfos;
	SamplingMutithread::runAlgMultiTask(tmp, solInfos, 0.5, 0.5, proname, params, algName, params, numRun);
	sols.clear();
	sols.resize(numRun);
	for (int idx(0); idx < tmp.size(); ++idx) {
		auto& cursol = tmp[idx];
		auto& curinfo = solInfos[idx];
		if (sols[curinfo->m_runId].size() <= curinfo->m_iter) {
			sols[curinfo->m_runId].resize(curinfo->m_iter + 1);
		}
		if (sols[curinfo->m_runId][curinfo->m_iter].size() <= curinfo->m_indiId) {
			sols[curinfo->m_runId][curinfo->m_iter].resize(curinfo->m_indiId + 1);
		}
		sols[curinfo->m_runId][curinfo->m_iter][curinfo->m_indiId] = cursol;
	}




	//sols = ofec::SamplingMutithread::runAlgMultiTask(0.5, 0.5, proname, params, algName, params, numRun);

}



void calculateNBN(
	ofec::nbn::NBNinfo& nbnInfo,
	//const std::string& saveDir,
	//const std::string& saveDir2,
	//const std::string& filename,
	ofec::Environment* env,
	ofec::Random* rnd
) {

	std::string proname = "TSP";
	auto pro = env->problem();
	//auto rnd = std::make_shared<ofec::Random>(0.5);
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
	nbnInfo.calculateNBNhnswEqualBetterRandom(env, rnd);
	//udpateNBNinfoIteration(nbnInfo,
	//	1e3, env, rnd);

	//nbnInfo.outputNBNinfo(saveDir, filename);
	//nbnInfo.outputVecSolutions<int>(saveDir, filename, env);


	//{

	//	//	auto filename = tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK);
	//	std::vector<int> solIds(nbnInfo.m_belong.size());
	//	for (int idx(0); idx < solIds.size(); ++idx) {
	//		solIds[idx] = idx;
	//	}
	//	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	//	ofec::transferNBNData(nbn_data, solIds, nbnInfo.m_belong, nbnInfo.m_dis2parent, nbnInfo.m_vFitness);

	//	ofec::TreeGraphSimple nbn_graph;
	//	nbn_graph.setNBNdata(nbn_data);
	//	//nbn_graph.modifyBestSols(rnd.get(), env->problem(), nbnInfo.solbases);
	//	nbn_graph.calNetwork(tsp_pro->numberVariables());


	//	std::string savepath = saveDir2 + filename + "_network.txt";
	//	std::ofstream out(savepath);
	//	nbn_graph.outputNBNnetwork(out);
	//	out.close();
	//	ouputNBNdata(saveDir2 + filename + "_nbn.txt", nbn_graph.get_nbn_data());
	//}



}



void calTask(const std::string& proReadDir,
	const std::string& lonRealDir,
	const std::string& saveDir,
	const std::string& tspfilename) {
	using namespace ofec;
	using namespace std;
	using namespace chrono;
	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));

	std::cout << "calculating filename\t" << tspname << std::endl;

	std::shared_ptr<Environment> env;
	genenrateTSPenv(proReadDir, tspfilename, env);
	int numSamples = 2e6;
	int nbnSamples = 1e6;
	numSamples = nbnSamples * 2;

	std::string proname = "TSP";
	auto pro = env->problem();
	auto rnd = std::make_shared<ofec::Random>(0.5);
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);

	auto eval_fun =
		[](ofec::SolutionBase& sol, ofec::Environment* env) {
		using namespace ofec;
		sol.evaluate(env, false);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};

	//std::shared_ptr<ofec::SolutionBase> bestSol;
	lon::LonInfo lon;
	lon::LonInfo filterLon;
	std::vector<std::vector<int>> lonSols;
	std::vector<int> lonSolNodeFit;
	std::vector<int> lonOriginId2newId;


	{

		std::string cursavePath = lonRealDir + tspname;
		cursavePath += "/";


		inputLon(cursavePath + tspname + "_origin_lon.txt", lon);
		lon::filterLON(lon, filterLon, 0.001, lonOriginId2newId);
		ofec::ParameterVariantStream variantStr;

		{
			std::stringstream buf;
			std::ifstream in(cursavePath + tspname + "_sol.txt");
			buf << in.rdbuf();
			in.close();
			ofec::variants_stream::stringstream2parameterStream(buf, variantStr);

		}

		size_t lonSolSize(0);
		variantStr >> lonSolSize;
		lonSols.resize(lonSolSize);
		for (auto& it : lonSols) {
			variantStr >> it;
		}
		//variantStr << lon_sampling.sols();
		variantStr >> lonSolNodeFit;



	}




	ofec::nbn::NBNinfo nbnInfo;
	for (auto& it : lonSols) {
		if (it.empty())continue;

		std::shared_ptr<ofec::SolutionBase> cursol(pro->createSolution());
		//cursol->initialize(env.get(), rnd);
		auto cursolx = it;
		for (auto& it : cursolx) --it;
		auto& cursolt = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*cursol);
		cursolt.variable().vect() = cursolx;
		//eval_fun(cursolt, env.get());
		nbnInfo.m_sols.push_back(cursol);
	}

	for (auto& it : nbnInfo.m_sols) {
		nbnInfo.m_solbases.push_back(it.get());
	}

	UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.m_solbases, env.get(), eval_fun);


	std::vector<int> bestSolIds(2, -1);
	std::vector<double> betSolFit(2, std::numeric_limits<double>::lowest());
	std::vector<int> bestFunnelId(2);
	std::vector<SolutionBase*>sols(2);
	for (int idx(0); idx < filterLon.m_node_funnel_idxs.size(); ++idx) {
		if (filterLon.m_node_funnel_idxs[idx] >= 0) {
			auto funId = filterLon.m_node_funnel_idxs[idx];
			auto solId = filterLon.m_nodeId[idx];
			double value = filterLon.m_node_fit[idx] + nbnInfo.m_solbases[solId]->fitness();
			//std::cout << "value\t" << value << std::endl;
			if (nbnInfo.m_solbases[solId]->fitness() > betSolFit[funId]) {
				bestSolIds[funId] = solId;
				betSolFit[funId] = nbnInfo.m_solbases[solId]->fitness();
				bestFunnelId[funId] = idx;
				sols[funId] = nbnInfo.m_solbases[solId];
			}
		}
	}

	{
		std::cout << "solDis\t" << sols.front()->variableDistance(sols.back()->variableBase(), env.get()) << std::endl;

	}

	{
		std::ofstream out(saveDir + tspname + "_lonBestFunnelIds.txt");
		ofec::matlab::outputVector(out, bestFunnelId);
		out.close();
	}




	calculateNBN(nbnInfo, env.get(), rnd.get());

	std::vector<int> showIdS = bestSolIds;

	//for (auto& it : bestSolIds) {
	//	if (nbnInfo.m_belong[it] != -1) {
	//		showIdS.push_back(nbnInfo.m_belong[it]);
	//	}
	////	std::cout << it << "\t";
	//}
	std::cout << std::endl;


	{
		ofec::nbn::NBNinfo curInfo;
		curInfo.subset(nbnInfo, showIdS);


		{
			outputNBNdata(std::cout, curInfo.m_solIds, curInfo.m_belong, curInfo.m_dis2parent, curInfo.m_vFitness);;
		}
	}



	{
		int updateSolId(0);
		if (sols.front()->fitness() > sols.back()->fitness()) {
			updateSolId = bestSolIds.back();
		}
		else {
			updateSolId = bestSolIds.front();
		}


		ofec::nbn::updateNearsetBetterSolution(updateSolId, nbnInfo.m_solbases, nbnInfo.m_belong, nbnInfo.m_dis2parent, nbnInfo.m_vFitness, env.get(), rnd.get());

	}


	{
		ofec::nbn::NBNinfo curInfo;
		curInfo.subset(nbnInfo, showIdS);

		std::cout << "update nbn NewInfo" << std::endl;
		{
			outputNBNdata(std::cout, curInfo.m_solIds, curInfo.m_belong, curInfo.m_dis2parent, curInfo.m_vFitness);;
		}
	}




	//for (auto& it : bestSolIds) {
	//	std::cout << nbnInfo.m_dis2parent[it] << "\t";
	//}
	//std::cout << std::endl;
	// 
	{
		std::vector<int> solIds;

		ofec::nbn::filterBestSolutions(solIds, nbnInfo, 0.9);
		for (auto& it : solIds) {
			std::cout << it << "\t";
		}
		std::cout << std::endl;
		//std::vector<ofec::SolutionBase*> centers;

	}
	//ofec::nbn::NBNinfo optInfo;

	//optInfo.subset(nbnInfo, solIds);

	//std::vector<int> centerIds;
	//{

	//	std::vector<int> solIds;
	//	ofec::nbn::getBestSolIds(optInfo, solIds);
	//	centerIds.push_back(optInfo.m_solIds[solIds.back()]);
	//	centers.push_back(optInfo.m_solbases[solIds.back()]);

	//	ofec::nbn::getSecondFurtherSolIds(optInfo, solIds);
	//	centerIds.push_back(optInfo.m_solIds[solIds.back()]);
	//	centers.push_back(optInfo.m_solbases[solIds.back()]);
	//	
	//}



	//bool flag_true(true);
	//for (int idx(0); idx < centers.size(); ++idx) {
	//	int newId = lonOriginId2newId[centerIds[idx]];
	//	int funnelId = filterLon.m_node_funnel_idxs[newId];
	//	if (funnelId >= 0 && betSolFit[funnelId] == centers[idx]->fitness()) {

	//	}
	//	else {
	//		flag_true = false;
	//		break;
	//	}
	//}

	//{
	//	std::ofstream out(saveDir + tspname + "_matlabBestFunnelIds.txt");
	//	ofec::matlab::outputVector(out, bestFunnelId);
	//	ofec::matlab::outputVector(out, centerIds);
	//	out.close();
	//}

	//if (flag_true) {
	//	std::cout << "yes, it is the further optima\t" << std::endl;
	//}
	//else {
	//	std::cout << "no, it is not the further optima\t" << std::endl;
	//}

	//{
	//	ofec::TreeGraphSimple nbn_graph;
	//	nbnInfo.trasferToGraph(nbn_graph, tsp_pro->numberVariables());
	//	std::string savepath = saveDir + tspname+ "_network.txt";
	//	std::ofstream out(savepath);
	//	nbn_graph.outputNBNnetwork(out);
	//	out.close();
	////	ouputNBNdata(savepath + tspname + "_nbn.txt", nbn_graph.get_nbn_data());
	//
	//}
	//



}






void runTask() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
	//	ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
		//ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";


	std::string readDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	readDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_data_eax_sampling_1000000_final/";
	//	saveDir = "F:/DiaoYiya/nbn_com_exp/nbn_sampling/eax_lkh_total_sortedX_2";
	//std::filesystem::create_directories(saveDir);

	readDir = "V:/DiaoYiya/code/data/ofec-data/paper_com_experiment_data/tsp_typical_comparison_exp_read/";
	readDir = "F:/DiaoYiya/nbn_com_exp_final/tsp_lon_data/";

	auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp3/nbn_lonData_funnelInfo/";
	std::filesystem::create_directories(saveDir2);

	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames =
	{ "5955.tsp",  "u574.tsp" , "2202.tsp", "1281.tsp" ,"6702.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", };
	filenames =
	{ "5955.tsp",  "u574.tsp" , "2202.tsp" };

	//std::reverse(filenames.begin(), filenames.end());

	for (auto& it : filenames) {
		if (it == "u574.tsp") {
			calTask(dir2, readDir, saveDir2, it);
		}
		else {
			calTask(dir1, readDir, saveDir2, it);
		}
	}





}



namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run(int argc, char* argv[]) {

		using namespace ofec;
		using namespace std;

		registerInstance();


		runTask();

	}


}