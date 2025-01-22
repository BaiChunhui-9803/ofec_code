
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
#include "../utility/nbn_visualization/nbn_fla/tsp_related/nbn_fla_tsp.h"
#include "../utility/nbn_visualization/instance/travelling_salemans_problem/nbn_modify_tsp.h"
#include "../utility/matlab/matlab_utility.h"

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




void calculateNBN(
	ofec::nbn::NBNinfo& nbnInfo,
	const std::string& saveDir,
	const std::string& saveDir2,
	const std::string& filename,
	ofec::Environment* env,
	ofec::Random* rnd
) {

	std::string proname = "TSP";
	auto pro = env->problem();
	//auto rnd = std::make_shared<ofec::Random>(0.5);
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
	nbnInfo.calculateNBNhnswEqualBetterRandom(env, rnd);
	nbnInfo.outputNBNinfo(saveDir, filename);
	nbnInfo.outputVecSolutions<int>(saveDir, filename, env);


	{

		nbnInfo.m_vFitness = UTILITY::valuesToSortedValues(nbnInfo.m_vFitness);
		//	auto filename = tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK);
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


void inputTrait(const std::string& saveDir, const std::string& filename, std::vector<std::vector<int>>& traits) {


	ofec::ParameterVariantStream paramsStream;
	ofec::variants_stream::inputFromFile(paramsStream, saveDir + filename + ".txt");
	size_t dataSize;
	paramsStream >> dataSize;
	traits.resize(dataSize);
	for (auto& it2 : traits) {
		paramsStream >> it2;
	}

}

void outputTrait(const std::string& saveDir, const std::string& filename, std::vector<std::vector<int>>& traits) {
	ofec::ParameterVariantStream paramsStream;
	paramsStream << traits.size();
	for (const auto& it2 : traits) {
		paramsStream << it2;

	}
	ofec::variants_stream::outputToFile(paramsStream, saveDir + filename + ".txt");
	{
		std::ofstream out(saveDir + filename + "_matlab.txt");
		auto& curTrait = traits;
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
		out.close();
	}

}

void filterToLDEE(
	const std::string& saveDri,
	const std::string& tspname,
	ofec::nbn::NBNinfo& nbnInfo, std::vector<ofec::SolutionBase*>& sols, std::vector<int>& toSolId) {
	using namespace ofec;
	using namespace std;
	//ofec::nbn::tsp::filterUniqueSols(sols, nbnInfo.m_solbases, toSolId, tsp_pro, rnd.get());


	{

		ofec::ParameterVariantStream paramsStream;
		paramsStream << toSolId;
		auto filepath = saveDri + tspname + "_mapSolId.txt";
		variants_stream::outputToFile(paramsStream, filepath);
	}

	{
		int sqrtSize = sqrt(nbnInfo.m_solbases.size());
		sqrtSize = sqrtSize * sqrtSize;
		std::cout << "origin size\t" << sols.size() << "\tcutsize\t" << sqrtSize << std::endl;
		nbnInfo.m_solbases.resize(sqrtSize);
		std::vector<std::vector<int>> tsp_sols;

		for (int idx(0); idx < nbnInfo.m_solbases.size(); ++idx) {
			auto& tspsol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*nbnInfo.m_solbases[idx]);
			tsp_sols.push_back(tspsol.variable().vect());
		}
		std::ofstream out(saveDri + tspname + "IDEE_solutions.txt");
		EAX_TSP_trait::outputIDEEsolution(out, tsp_sols);
		out.close();
	}
}


ofec::nbn::NBNinfo nbnInfo;
void readEAX_LKH(const std::string& proDir,
	const std::string& eaxlkhDir,
	const std::string& lonDir,
	const std::string& saveDri,
	const std::string& tspfilename) {
	using namespace ofec;
	using namespace std;
	using namespace chrono;
	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));
	std::cout << "calculating filename\t" << tspname << std::endl;
	std::shared_ptr<Environment> env;
	genenrateTSPenv(proDir, tspfilename, env);

	std::string proname = "TSP";
	auto pro = env->problem();
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
	auto rnd = std::make_shared<ofec::Random>(0.5);

	ofec::nbn::NBNinfo eaxInfo;
	{
		std::cout << "inputNBNinfo solutions of eax lkh by multithread" << std::endl;
		//eaxInfo.inputNBNinfo(eaxlkhDir, tspname + "_eax_lkh");
		eaxInfo.inputVecSolutions<int>(eaxlkhDir, tspname, env.get());
	}

	ofec::nbn::NBNinfo lonNBNinfo;
	{
		ofec::nbn::tsp::TspSolutionsToNBNinfo(lonNBNinfo, lonDir, tspname + "_lonSamples", env.get());
		//lonNBNinfo.tsp(lonDir, tspname + "_lonSamples", env.get());
	}

	//int lonSolSize = lonNBNinfo.solbases.size();



	std::vector<ofec::SolutionBase*>sols;
	for (auto& it : lonNBNinfo.m_solbases) {
		sols.push_back(it);
	}
	for (auto& it : eaxInfo.m_solbases) {
		sols.push_back(it);
	}
	auto eval_fun =
		[](ofec::SolutionBase& sol, ofec::Environment* env) {
		using namespace ofec;
		sol.evaluate(env, false);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};

	UTILITY::evaluateRandomSolutionsMultiThreads(sols, env.get(), eval_fun);



	std::vector<int> toSolId;
	ofec::nbn::tsp::filterUniqueSols(sols, nbnInfo.m_solbases, toSolId, tsp_pro, rnd.get());

	//{
	//	auto filepath = saveDri + tspname + "_mapSolId.txt";
	//	ofec::ParameterVariantStream paramsStream;
	//	variants_stream::inputFromFile(paramsStream, filepath);
	//	paramsStream >> toSolId;


	//}

	{
		//int maxSolId(0);
		//ofec::calMax(toSolId,maxSolId);
		//nbnInfo.m_solbases.resize(maxSolId + 1);
		//for (int idx(0); idx < nbnInfo.m_solbases.size(); ++idx) {
		//	nbnInfo.m_solbases[idx] = sols[toSolId[idx]];
		//}


		int sqrtSize = sqrt(nbnInfo.m_solbases.size());
		sqrtSize = sqrtSize * sqrtSize;
		std::cout << "origin size\t" << sols.size() << "\tcutsize\t" << sqrtSize << std::endl;
		nbnInfo.m_solbases.resize(sqrtSize);

		//int sqrtSolSize = sqrt(toSolId.size());
		//sqrtSolSize = sqrtSolSize * sqrtSolSize;
		//std::cout << "origin size\t" << sols.size() << "\tcutsize\t" << sqrtSolSize << std::endl;
		//nbnInfo.m_solbases.resize(sqrtSolSize);



		{
			std::vector<double> fitness;
			int solId(0);
			for (auto& it : nbnInfo.m_solbases) {

				//	 std::cout << "solId\t" << solId++ << std::endl;
				if (it == nullptr) std::cout << "empty ptr\t" << solId << std::endl;
				fitness.push_back(it->fitness());
				++solId;
			}


			double maxFit(0);
			calMax(fitness, maxFit);
			std::vector<int> solIds;
			for (int idx(0); idx < fitness.size(); ++idx) {
				if (fitness[idx] == maxFit) {
					solIds.push_back(idx);
				}
			}

			{
				std::ofstream out(saveDri + tspname + "_bestSolIds.txt");
				ofec::matlab::outputVector(out, solIds);
				out.close();
			}

			ofec::ParameterVariantStream paramsStream;
			paramsStream << fitness;
			ofec::variants_stream::outputToFile(paramsStream, saveDri + tspname + "_fitness.txt");

			{
				fitness = UTILITY::valuesToSortedValues(fitness);

				std::ofstream out(saveDri + tspname + "_fitness_matlab.txt");
				ofec::matlab::outputVector(out, fitness);
				out.close();

			}


			//int runId(0);

			for (int runId(0); runId < 30; ++runId) {
				std::vector<std::vector<int>> trait;
				inputTrait(eaxlkhDir, tspname + "_lkh_trait_" + std::to_string(runId), trait);
				for (auto& it : trait) {
					for (auto& it2 : it) {
						it2 = toSolId[it2 + lonNBNinfo.m_solbases.size()];
					}
				}

				outputTrait(saveDri, tspname + "_lkh_trait_" + std::to_string(runId), trait);
			}


			for (int runId(0); runId < 30; ++runId) {
				std::vector<std::vector<int>> trait;
				inputTrait(eaxlkhDir, tspname + "_eax_trait_" + std::to_string(runId), trait);
				for (auto& it : trait) {
					for (auto& it2 : it) {
						it2 = toSolId[it2 + lonNBNinfo.m_solbases.size()];
					}
				}

				outputTrait(saveDri, tspname + "_eax_trait_" + std::to_string(runId), trait);
			}


		}



	}


}


void runTask() {
	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//\\172.24.34.11\share\2018\diaoyiya\ofec-data\paper_com_experiment_data\totalTsp\total_data_idee
	ofec::g_working_directory = "//172.24.34.11/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";



//	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	auto eax_dir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_more_hnsw_rnd_3_linux/";

	auto lonDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_data_lonSample/";

	auto ldeeSave = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/total_data_idee/";
	std::filesystem::create_directories(ldeeSave);
	//std::filesystem::create_directories(saveDir2);

	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";


	std::vector<std::string> filenames =
	{ "2202.tsp","1281.tsp" ,"u574.tsp" ,"6717.tsp" , "6702.tsp" , "7310.tsp", "9225.tsp", "5955.tsp", };
	//	filenames.clear();
	//	filenames = { "2202.tsp" };

		//	std::reverse(filenames.begin(), filenames.end());


	std::cout << "IDEE solutions output2" << std::endl;
	filenames =
	{ "5955.tsp",  "u574.tsp" , "2202.tsp" };


	//filenames =
	//{ "2202.tsp" };
	for (auto& it : filenames) {
		if (it == "u574.tsp") {
			readEAX_LKH(dir2, eax_dir, lonDir, ldeeSave, it);
		}
		else {
			readEAX_LKH(dir1, eax_dir, lonDir, ldeeSave, it);
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