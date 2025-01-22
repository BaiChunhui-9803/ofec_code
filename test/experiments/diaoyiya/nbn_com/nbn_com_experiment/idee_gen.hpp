
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

#include "../utility/nbn_visualization/nbn_fla/nbn_fla.h"
#include "../instance/algorithm/combination/eax_tsp/eax_tsp_basin/eax_tsp_basin_multipop.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
#include "../instance/algorithm/combination/LKH_origin/INCLUDE/LKH.h"


#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>


#include "../utility/visualization/idee_trait.h"



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



void filterUniqueSols(const std::vector<ofec::SolutionBase*>& originSol,
	std::vector<ofec::SolutionBase*>& filteredSols,
	std::vector<int>& toSolId,
	ofec::TravellingSalesman* tsp_pro, ofec::Random* rnd
) {

	ofec::TravellingSalesman::HashSolutionMap solMap;
	solMap.initialize(rnd, tsp_pro->numberVariables() + 10);
	toSolId.resize(originSol.size());
	filteredSols.clear();
	for (int idx(0); idx < originSol.size(); ++idx) {
		toSolId[idx] = solMap.getSolId(*originSol[idx]);
		if (toSolId[idx] == filteredSols.size()) {
			filteredSols.push_back(originSol[idx]);
		}
	}

}

void outputToFile(ofec::ParameterVariantStream& paramsStream, const std::string& filepath) {
	std::stringstream buf;
	ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
	std::ofstream out(filepath);
	out << buf.rdbuf();
	out.close();
}


void calTask(const std::string& readDir,
	const std::string& saveDir,
	const std::string& tspfilename) {
	using namespace ofec;
	using namespace std;
	using namespace chrono;
	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));

	std::cout << "calculating filename\t" << tspname << std::endl;

	std::shared_ptr<Environment> env;
	genenrateTSPenv(readDir, tspfilename, env);
	int numSamples = 2e6;
	int nbnSamples = 2e5;
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


	{


		std::cout << "generate solutions by multithread EAX-TSP-sampling" << std::endl;
		std::string algName = "EAX-TSP-sampling";
		ofec::ParameterMap params;
		//params["problem name"] = std::string("TSP");
		params["dataFile1"] = tspfilename;
		params["dataDirectory"] = readDir;
		//params["algorithm name"] = std::string("EAX-TSP-sampling");
		auto sols = ofec::SamplingMutithread::runAlgMultiTask(0.5, 0.5, proname, params, algName, params, 30);
		// auto sols = ofec::SamplingMutithread::runAlgMultiTask(pro.get(), algName, params, numRun);

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

		std::vector<int> filterSolId;
		std::vector<ofec::SolutionBase*> totalSols;
		filterUniqueSols(solbases, totalSols, filterSolId, tsp_pro, rnd.get());

		UTILITY::evaluateRandomSolutionsMultiThreads(totalSols, env.get(), eval_fun);

		for (auto& it : solIds) {
			for (auto& it2 : it) {
				for (auto& it3 : it2) {
					it3 = filterSolId[it3];
				}
			}
		}
		std::vector<ofec::Real> solFit;
		for (auto& it : totalSols) {
			solFit.push_back(it->fitness());
		}

		{
			ofec::ParameterVariantStream paramsStream;
			paramsStream << solIds.size();
			for (const auto& it : solIds) {
				paramsStream << it.size();
				for (const auto& it2 : it) {
					paramsStream << it2;
				}
			}

			paramsStream << solFit;
			outputToFile(paramsStream, saveDir + tspname + "_eaxRun_solIds_fitness.txt");

		}


		ofec::EAX_TSP_trait trait;
		for (int idx(0); idx < totalSols.size(); ++idx) {
			auto& tspsol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*totalSols[idx]);
			trait.m_sols.push_back(tspsol.variable().vect());
		}

		trait.outputIDEESols(saveDir + tspname + "IDEE_solutions.txt");
	}




}






void runTask() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//		ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//		ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
		//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
	//		ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/NBN_datacode_total/data";

	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/idee_datas/";

	std::filesystem::create_directories(saveDir);


	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames =
	{ "2202.tsp","1281.tsp" ,"6702.tsp" ,"u574.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", "5955.tsp", };


	//	std::reverse(filenames.begin(), filenames.end());

	for (auto& it : filenames) {
		if (it == "u574.tsp") {
			calTask(dir2, saveDir, it);
		}
		else {
			calTask(dir1, saveDir, it);
		}
	}

}



namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run() {

		using namespace ofec;
		using namespace std;

		registerInstance();


		runTask();

	}


}