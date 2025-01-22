
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
#include "../utility/nbn_visualization/nbn_fla/tsp_related/nbn_fla_tsp.h"\


#include "../utility/nbn_visualization/nbn_calculator/nbn_improve_accuracy_multithread.h"
#include "../utility/nbn_visualization/nbn_calculator/nbn_modifiy_equal.h"







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



void calTaskOneTask(
	const std::string& saveDir,
	const std::string& saveDir2,
	const std::string& tspfilename,
	int neighborK,
	int nbnSamples,
	ofec::Environment* env,
	ofec::Random* rnd) {


	using namespace ofec;
	using namespace std;
	using namespace nbn;
	using namespace tsp;


	int numSamples = 2e6;
	//int nbnSamples = 1e6;
	numSamples = nbnSamples * 2;
	auto eval_fun =
		[](ofec::SolutionBase& sol, ofec::Environment* env) {
		using namespace ofec;
		sol.evaluate(env, false);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};
	auto pro = env->problem();
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);

	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));
	std::string filename;
	std::vector<std::shared_ptr<ofec::SolutionBase>> samples;

	if (neighborK < 0) {
		filename = tspname + "_randomSamples";
	}
	else {
		filename = tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK);
	}

	std::string filepath = { nbn::NBNinfo::getSolutionsFileName(saveDir,filename) };
	if (std::filesystem::exists(filepath))
	{
		ofec::nbn::NBNinfo nbnInfo;
		nbnInfo.inputVecSolutions<int>(saveDir, filename, env);
		nbnInfo.inputNBNinfo(saveDir, filename);

		ofec::nbn::modifyEqualSolutions(nbnInfo, env, rnd);

		{
			std::cout << "outputNBNinfo solutions by multithread" << std::endl;
			nbnInfo.outputNBNinfo(saveDir2, filename);
		}


		{

			//auto filename = tspname + "_randomSamples";
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




void calTaskTotal(const std::string& readDir,
	const std::string& saveDir,
	const std::string& saveDir2,
	const std::string& tspfilename,
	const std::vector<int>& tasks
) {
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

	for (auto& curTask : tasks) {
		calTaskOneTask(
			saveDir, saveDir2, tspfilename, curTask, nbnSamples, env.get(), rnd.get());
	}


}





void runTask2() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
	ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";

	//	ofec::g_working_directory = "//172.24.24.151/e/DiaoYiya/code/data/ofec-data";

	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp2/nbn_data_eax_lkh_compare/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread/";
	std::filesystem::create_directories(saveDir);

	auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread_equalBetter/";
	std::filesystem::create_directories(saveDir2);

	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames /*=
	{ "2202.tsp", "u574.tsp" ,"5955.tsp", "1281.tsp" ,"6702.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", }*/;

	filenames = { "u574.tsp",  "5955.tsp" ,   "2202.tsp", };

	std::vector<std::pair<std::string, std::vector<int>>> tasks;


	//{3,6,12,}
	std::vector<int> neighbors = { -1, 100, 50,25,12,6,3 };

	tasks.push_back({ "5955.tsp" , std::vector<int>({ -1, 100, 50,25,12,6,3}) });
	tasks.push_back({ "u574.tsp" , std::vector<int>({ /*-1,*/ 100, 50,25/*,12*/,6/*,3*/}) });

	//	tasks.push_back({ "2202.tsp" , std::vector<int>({ -1, 100, 50,25,12,6,3}) });

	tasks.push_back({ "u574.tsp" , std::vector<int>({ -1,12,3}) });


	std::reverse(tasks.begin(), tasks.end());

	for (auto& it : tasks) {
		//	it.second = neighbors;

		std::reverse(it.second.begin(), it.second.end());
	}

	std::cout << "hnsw single thread " << std::endl;

	for (auto& it : tasks) {
		if (it.first == "u574.tsp") {
			calTaskTotal(dir2, saveDir, saveDir2, it.first, it.second);
			//calTask(dir2, saveDir, saveDir2, it);
		}
		else {
			calTaskTotal(dir1, saveDir, saveDir2, it.first, it.second);
			//calTask(dir1, saveDir, saveDir2, it);
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


		runTask2();

	}


}