
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
#include "../utility/nbn_visualization/calculate_algorithm/nbn_grid_tree_division.h"
#include "../utility/nbn_visualization/calculate_algorithm/nbn_grid_tree_division_allSize.h"
#include  "../utility/function/custom_function.h"
#include "../utility/function/general_function.h"
#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator.h"
//#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include "interface.h"
//#include "../utility/function/custom_function.h"



#include "../core/algorithm/population.h"
//#include "../utility/nondominated_sorting/filter_sort.h"
#include "../instance/algorithm/visualize/sampling/instance/tsp/sampling_eax_tsp.h"
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





void calTask(const std::string& readDir,
	const std::string& saveDir,
	const std::string& saveDir2,


	const std::string& tspfilename) {
	using namespace ofec;
	using namespace std;
	using namespace nbn;
	using namespace tsp;

	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));

	std::cout << "calculating filename\t" << tspname << std::endl;

	std::shared_ptr<Environment> env;
	genenrateTSPenv(readDir, tspfilename, env);
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

	std::shared_ptr<ofec::SolutionBase> bestSol;



	//{

	//	std::string filepath = { nbn::NBNinfo::getSolutionsFileName(saveDir, tspname + "_randomSamples") };
	//	bool filepathExists = std::filesystem::exists(filepath);
	//	if (!filepathExists) {
	//		std::vector<std::shared_ptr<ofec::SolutionBase>> samples;
	//		std::cout << "generate solutions by multithread randomSample" << std::endl;
	//		ofec::nbn::tsp::sampleSolsRandomMultiThread(samples, numSamples, env.get(), rnd.get());
	//		std::vector<SolutionBase*> solbases;
	//		for (auto& it : samples) {
	//			solbases.push_back(it.get());
	//		}

	//		NBNinfo nbnInfo;
	//		std::cout << "filterUniqueSols solutions by multithread" << std::endl;
	//		{
	//			std::vector<int> filterSolId;
	//			filterUniqueSols(solbases, nbnInfo.m_solbases, filterSolId, tsp_pro, rnd.get());
	//		}
	//		nbn::tsp::sortEachSolutionX(nbnInfo.m_solbases);
	//		UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.m_solbases, env.get(), eval_fun);

	//		if (nbnInfo.m_solbases.size() > nbnSamples)
	//			nbnInfo.m_solbases.resize(nbnSamples);
	//		std::cout << "calculateNBN solutions by multithread" << std::endl;
	//		nbnInfo.calculateNBNhnswEqualBetter(env.get(), rnd.get());

	//		std::cout << "outputNBNinfo solutions by multithread" << std::endl;
	//		nbnInfo.outputNBNinfo(saveDir, tspname + "_randomSamples");
	//		nbnInfo.outputVecSolutions<int>(saveDir, tspname + "_randomSamples", env.get());


	//		{

	//			auto filename = tspname + "_randomSamples";
	//			std::vector<int> solIds(nbnInfo.m_belong.size());
	//			for (int idx(0); idx < solIds.size(); ++idx) {
	//				solIds[idx] = idx;
	//			}
	//			std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	//			ofec::transferNBNData(nbn_data, solIds, nbnInfo.m_belong, nbnInfo.m_dis2parent, nbnInfo.m_vFitness);

	//			ofec::TreeGraphSimple nbn_graph;
	//			nbn_graph.setNBNdata(nbn_data);
	//			//nbn_graph.modifyBestSols(rnd.get(), env->problem(), nbnInfo.solbases);
	//			nbn_graph.calNetwork(tsp_pro->numberVariables());


	//			std::string savepath = saveDir2 + filename + "_network.txt";
	//			std::ofstream out(savepath);
	//			nbn_graph.outputNBNnetwork(out);
	//			out.close();
	//			ouputNBNdata(saveDir2 + filename + "_nbn.txt", nbn_graph.get_nbn_data());
	//		}
	//	}
	//}


	{


		std::cout << "generate solutions by multithread EAX-TSP-sampling" << std::endl;
		std::string algName = "EAX-TSP-sampling";
		ofec::ParameterMap params;
		//params["problem name"] = std::string("TSP");
		params["dataFile1"] = tspfilename;
		params["dataDirectory"] = readDir;
		//params["algorithm name"] = std::string("EAX-TSP-sampling");
		//auto sols = ofec::SamplingMutithread::runAlgMultiTask(0.5, 0.5, proname, params, algName, params, 30);
		// auto sols = ofec::SamplingMutithread::runAlgMultiTask(pro.get(), algName, params, numRun);


		std::vector<std::vector<std::vector<std::shared_ptr<ofec::SolutionBase>>>> sols;
		runEax(readDir, tspfilename, sols);

		std::vector<ofec::SolutionBase*> solbases;
		std::vector<std::vector<std::vector<int>>> solIds(sols.size());

		for (int idx(0); idx < sols.size(); ++idx) {
			solIds[idx].resize(sols[idx].size());
			for (int idy(0); idy < sols[idx].size(); ++idy) {
				solIds[idx][idy].resize(sols[idx][idy].size());
				for (int idz(0); idz < sols[idx][idy].size(); ++idz) {
					solIds[idx][idy][idz] = solbases.size();
					solbases.push_back(sols[idx][idy][idz].get());
				}
			}
		}


		NBNinfo nbnInfo;

		{
			std::vector<int> filterSolId;
			filterUniqueSols(solbases, nbnInfo.m_solbases, filterSolId, tsp_pro, rnd.get());
		}
		nbn::tsp::sortEachSolutionX(nbnInfo.m_solbases);
		UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.m_solbases, env.get(), eval_fun);

		bestSol = nbnInfo.m_solbases.front()->getSharedPtr();
		for (auto& it : nbnInfo.m_solbases) {
			if (bestSol->fitness() < it->fitness()) {
				bestSol = it->getSharedPtr();
			}
		}



	}


	std::vector<int> neighbors = { 8 };
	//std::reverse(neighbors.begin(), neighbors.end());

	auto& bestTspSol = dynamic_cast<TravellingSalesman::SolutionType&>(*bestSol);
	for (auto& neighborK : neighbors) {

		std::cout << "generate solutions by multithread sampling neighborK\t" << neighborK << std::endl;


		std::string filepath = { NBNinfo::getSolutionsFileName(saveDir,tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK)) };
		bool filepathExists = std::filesystem::exists(filepath);

		if (!filepathExists) {
			std::vector<std::shared_ptr<ofec::SolutionBase>> samples;
			nbn::tsp::sampleSolsAroundTSPMultiThread_kSwap(bestTspSol, samples, neighborK, numSamples, env.get(), rnd.get());
			std::vector<SolutionBase*> solbases;
			for (auto& it : samples) {
				solbases.push_back(it.get());
			}

			NBNinfo nbnInfo;
			{
				std::vector<int> filterSolId;
				filterUniqueSols(solbases, nbnInfo.m_solbases, filterSolId, tsp_pro, rnd.get());
				//nbn::tsp::sortEachSolutionX(nbnInfo.m_solbases);
			}

			if (nbnInfo.m_solbases.size() > nbnSamples)
				nbnInfo.m_solbases.resize(nbnSamples);

			nbn::tsp::sortEachSolutionX(nbnInfo.m_solbases);

			UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.m_solbases, env.get(), eval_fun);
			nbnInfo.calculateNBNhnswEqualBetterRandom(env.get(), rnd.get());

			nbnInfo.outputNBNinfo(saveDir, tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK));
			nbnInfo.outputVecSolutions<int>(saveDir, tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK), env.get());


			{

				auto filename = tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK);
				std::vector<int> solIds(nbnInfo.m_belong.size());
				for (int idx(0); idx < solIds.size(); ++idx) {
					solIds[idx] = idx;
				}
				std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
				ofec::transferNBNData(nbn_data, solIds, nbnInfo.m_belong, nbnInfo.m_dis2parent, nbnInfo.m_vFitness);

				ofec::TreeGraphSimple nbn_graph;
				nbn_graph.setNBNdata(nbn_data);
				//nbn_graph.modifyBestSols(rnd.get(), env->problem(), nbnInfo.solbases);
				nbn_graph.calNetwork(tsp_pro->numberVariables());


				std::string savepath = saveDir2 + filename + "_network.txt";
				std::ofstream out(savepath);
				nbn_graph.outputNBNnetwork(out);
				out.close();
				ouputNBNdata(saveDir2 + filename + "_nbn.txt", nbn_graph.get_nbn_data());
			}



		}
	}

}




void runTask() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";

	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp2/nbn_data_eax_lkh_compare/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread_2/";
	std::filesystem::create_directories(saveDir);

	auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread_2_network/";
	std::filesystem::create_directories(saveDir2);

	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames /*=
	{ "2202.tsp", "u574.tsp" ,"5955.tsp", "1281.tsp" ,"6702.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", }*/;

	filenames = { "u574.tsp",  "5955.tsp" ,   "2202.tsp", };


	//std::reverse(filenames.begin(), filenames.end());


	std::cout << "calculating nbn subresion all task ver4" << std::endl;
	//readTspInstance(filenames);

	//size_t K = 7; // 假设我们想要删除前3个元素

	//if (K <= filenames.size()) {
	//	filenames.erase(filenames.begin(), std::next(filenames.begin(), K));
	//}
	// 
	// 
//	filenames.resize(31);
//	std::reverse(filenames.begin(), filenames.end());
	//std::reverse(filenames.begin(), filenames.end());

	for (auto& it : filenames) {
		if (it == "u574.tsp") {
			calTask(dir2, saveDir, saveDir2, it);
		}
		else {
			calTask(dir1, saveDir, saveDir2, it);
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