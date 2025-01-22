
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
#include "../utility/nbn_visualization/visualization/tree_graph_simple.h"
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
#include "../instance/problem/combination/travelling_salesman/tsp_offline_data/travelling_salesman_offline_data.h"

#include "../utility/nbn_visualization/nbn_fla/nbn_fla.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
#include "../instance/algorithm/combination/LKH_origin/INCLUDE/LKH.h"
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>

#include "../utility/idee/idee_trait.h"
#include "../utility/nbn_visualization/instance/travelling_salemans_problem/nbn_modify_tsp.h"

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




void readEAX_LKH(const std::string& proDir,
	const std::string& eaxlkhDir,
	const std::string& lonDir,
	const std::string& saveDri,
	const std::string& tspfilename) {
	using namespace ofec;
	using namespace std;
	using namespace nbn;

	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));
	std::cout << "calculating filename\t" << tspname << std::endl;
	std::shared_ptr<Environment> env;
	genenrateTSPenv(proDir, tspfilename, env);

	std::string proname = "TSP";
	auto pro = env->problem();
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
	auto rnd = std::make_shared<ofec::Random>(0.5);

	NBNinfo eaxInfo;
	{
		std::cout << "inputNBNinfo solutions by multithread" << std::endl;
		//eaxInfo.inputNBNinfo(eaxlkhDir, tspname + "_eax_lkh");

		eaxInfo.inputVecSolutions<int>(eaxlkhDir, tspname + "_eax_lkh", env.get());
	}

	std::vector<std::string> filenames = { "lkh_fail", "lkh_success", "eax_fail", "eax_success" };
	std::vector<std::vector<std::vector<int>>> m_traits(4);
	std::vector<std::string> task_names = { "lkh", "eax" };
	std::vector<int> successRate;

	{
		ofec::ParameterVariantStream paramsStream;
		ofec::variants_stream::inputFromFile(paramsStream, eaxlkhDir + tspname + "_trait_sucessRate.txt");
		size_t inputSize(0);
		paramsStream >> inputSize;
		m_traits.resize(inputSize);
		for (auto& it : m_traits) {
			paramsStream >> inputSize;
			it.resize(inputSize);
			for (auto& it2 : it) {
				paramsStream >> it2;
			}
		}
		paramsStream >> successRate;


	}




	NBNinfo lonNBNinfo;
	{

		tsp::TspSolutionsToNBNinfo(lonNBNinfo, lonDir, tspname + "_lonSamples", env.get());
		//	lonNBNinfo.inputTspSolutions(lonDir, tspname + "_lonSamples", env.get());
	}

	int lonSolSize = lonNBNinfo.m_solbases.size();
	NBNinfo nbnInfo;


	std::vector<ofec::SolutionBase*>sols;
	for (auto& it : lonNBNinfo.m_solbases) {
		sols.push_back(it);
	}

	for (auto& it : eaxInfo.m_solbases) {
		sols.push_back(it);
	}
	std::vector<int> toSolId;
	tsp::filterUniqueSols(sols, nbnInfo.m_solbases, toSolId, tsp_pro, rnd.get());



	auto eval_fun =
		[](ofec::SolutionBase& sol, ofec::Environment* env) {
		using namespace ofec;
		sol.evaluate(env, false);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};

	UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.m_solbases, env.get(), eval_fun);


	nbnInfo.calculateNBNhnswEqualBetterRandom(env.get(), rnd.get());
	{
		std::cout << "outputNBNinfo solutions by multithread" << std::endl;
		nbnInfo.outputNBNinfo(saveDri, tspname + "_lon_eax");
		nbnInfo.outputVecSolutions<int>(saveDri, tspname + "_lon_eax", env.get());
	}

	std::vector<int> solIds(nbnInfo.m_solbases.size());
	for (int idx(0); idx < solIds.size(); ++idx) {
		solIds[idx] = idx;
	}

	{
		nbnInfo.m_vFitness = UTILITY::valuesToSortedValues(nbnInfo.m_vFitness);
	}

	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	ofec::transferNBNData(nbn_data, solIds, nbnInfo.m_belong, nbnInfo.m_dis2parent, nbnInfo.m_vFitness);

	ofec::TreeGraphSimple nbn_graph;
	nbn_graph.setNBNdata(nbn_data);
	nbn_graph.modifyBestSols(rnd.get(), env->problem(), nbnInfo.m_solbases);
	nbn_graph.calNetwork(tsp_pro->numberVariables());


	std::string savepath = saveDri + tspname + "_network.txt";
	std::ofstream out(savepath);
	nbn_graph.outputNBNnetwork(out);
	out.close();
	ouputNBNdata(saveDri + tspname + "_nbn.txt", nbn_graph.get_nbn_data());


	{
		for (auto& it : m_traits) {
			for (auto& it2 : it) {
				for (auto& it3 : it2) {
					it3 = toSolId[it3 + lonSolSize];
				}
			}
		}
	}



	{
		std::ofstream out(saveDri + tspname + "_trait_sucessRate.txt");
		for (auto& it : task_names) {
			out << it << "\t";
		}
		out << std::endl;
		for (auto& it : successRate) {
			out << it << "\t";
		}
		out << std::endl;
		out.close();
	}



	{
		std::ofstream out(saveDri + tspname + "_traits.txt");
		for (auto& curTrait : m_traits) {
			std::map<int, int> mii;
			//	int maxIter(0);
			for (int itIter(0); itIter < curTrait.size(); ++itIter) {
				for (auto& it2 : curTrait[itIter]) {
					if (it2 > 0) {
						mii[it2] = std::max(mii[it2], itIter);
					}
				}
			}


			out << 2 << "\t" << mii.size() << std::endl;
			for (auto& it : mii) {
				out << it.first << "\t" << double(it.second) / double(curTrait.size()) << std::endl;
			}

		}

		out.close();
	}


	{
		double maxFit(0);
		auto& fitness = nbnInfo.m_vFitness;
		ofec::calMax(fitness, maxFit);
		std::vector<int> bestIds;
		for (int idx(0); idx < fitness.size(); ++idx) {
			if (fitness[idx] == maxFit) {
				std::cout << idx << "\t" << maxFit << std::endl;
				std::cout << "maxFit\t" << maxFit << std::endl;
				bestIds.push_back(idx);
			}
		}
		{
			std::ofstream out2(saveDri + tspname + "_bestIds.txt");
			out2 << 1 << "\t" << bestIds.size() << std::endl;
			for (auto& it : bestIds) {
				out2 << it << std::endl;
			}
			out2.close();
		}
	}

	//{
	//	std::ofstream out(saveDri + tspname + "_nodeFit.txt");
	//	auto& fitness = nbnInfo.vFitness;
	//	ofec::dataNormalize(fitness);
	//	out << 1 << "\t" << fitness.size() << std::endl;
	//	for (auto& it : fitness) {
	//		out << it << std::endl;
	//	}
	//	out.close();
	//}



}


void runTask() {
	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";



	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp3/eax_lkh_nbn_hnswRnd/";

	auto saveDir1 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp3/nbn_data_lonSample/";


	auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp3/eax_lkh_lon_nbndata_rnd/";
	std::filesystem::create_directories(saveDir);
	std::filesystem::create_directories(saveDir2);

	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";


	std::vector<std::string> filenames =
	{ "2202.tsp","1281.tsp" ,"u574.tsp" ,"6717.tsp" , "6702.tsp" , "7310.tsp", "9225.tsp", "5955.tsp", };
	//	filenames.clear();
	//	filenames = { "2202.tsp" };

		//	std::reverse(filenames.begin(), filenames.end());
	filenames =
	{ "5955.tsp",  "u574.tsp" , "2202.tsp" };
	for (auto& it : filenames) {
		if (it == "u574.tsp") {
			readEAX_LKH(dir2, saveDir, saveDir1, saveDir2, it);
		}
		else {
			readEAX_LKH(dir1, saveDir, saveDir1, saveDir2, it);
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