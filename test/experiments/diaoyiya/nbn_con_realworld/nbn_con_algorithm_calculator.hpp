
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
#include "../utility/function/custom_function.h"
#include "../utility/function/general_function.h"
#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator.h"
#include "interface.h"
#include <queue>
#include "../core/definition.h"
#include <deque>
#include <chrono>   
#include "../utility/hnsw/hnsw_nbn/neighbor_info.h"
#include "../utility/hnsw/hnsw_nbn/select_policy.h"
#include "../utility/hnsw/hnsw_nbn/node.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_model.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_nbn.h"
#include "../core/environment/environment.h"
#include "../core/parameter/variants_to_stream.h"


#include <filesystem>
#include <queue>
#include <set>

#include "../core/parameter/variants_to_stream.h"
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>



#include "../utility/nbn_visualization/instance/travelling_salemans_problem/nbn_modify_tsp.h"
#include "../utility/nbn_visualization/calculate_algorithm/nbn_iteration_multithread.h"
#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"


#include "../utility/hnsw/hnsw_nbn/neighbor_info.h"
#include "../utility/hnsw/hnsw_nbn/select_policy.h"
#include "../utility/hnsw/hnsw_nbn/node.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_model.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_nbn.h"
#include "../core/environment/environment.h"
#include "../core/parameter/variants_to_stream.h"
#include "../utility/function/custom_function.h"
#include "../core/parameter/variants_to_stream.h"

#include "../instance/algorithm/continuous/single_objective/global/cma_es/cma_es.h"

#include "../instance/algorithm/visualize/sampling_datatum/sampling_algorithm_multipop.h"
#include "../instance/algorithm/visualize/sampling_datatum/sampling_algorithm_multipop_multithread.hpp"
#include "../instance/algorithm/visualize/sampling_datatum/instance/continous/global/cmaes_sampling.h"


void calculateNBN(
	ofec::nbn::NBNinfo& nbnInfo,
	const std::string& saveDir,
	const std::string& saveDir2,
	const std::string& filename,
	ofec::Environment* env,
	ofec::Random* rnd
) {

	//std::string proname = "TSP";
	auto pro = env->problem();
	//auto rnd = std::make_shared<ofec::Random>(0.5);
	//auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);

	nbnInfo.calculateNBNhnswEqualBetterRandom(env, rnd);



	//udpateNBNinfoIteration(nbnInfo,
	//	1e3, env, rnd);
	std::cout << "outputNBNinfo" << std::endl;
	nbnInfo.outputNBNinfo(saveDir, filename);
	//	nbnInfo.outputVecSolutions<int>(saveDir, filename, env);

	{
		std::cout << "output solutions" << std::endl;
		ofec::ParameterVariantStream paramsStream;
		paramsStream << nbnInfo.m_solbases.size();
		for (auto& it : nbnInfo.m_solbases) {
			pro->solutionToParameterVariants(*it, paramsStream);
		}
		ofec::variants_stream::outputToFile(paramsStream, saveDir + filename + ".txt");
	}

	{


		nbnInfo.m_vFitness = UTILITY::valuesToSortedValues(nbnInfo.m_vFitness);

		//	auto filename = tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK);
		std::vector<int> solIds(nbnInfo.m_belong.size());
		for (int idx(0); idx < solIds.size(); ++idx) {
			solIds[idx] = idx;
		}
		std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
		ofec::transferNBNData(nbn_data, solIds, nbnInfo.m_belong, nbnInfo.m_dis2parent, nbnInfo.m_vFitness);
		std::cout << "setNBNdata" << std::endl;
		ofec::TreeGraphSimple nbn_graph;
		nbn_graph.setNBNdata(nbn_data);

		std::cout << "calNetwork" << std::endl;
		//nbn_graph.modifyBestSols(rnd.get(), env->problem(), nbnInfo.solbases);
		nbn_graph.calNetwork();


		std::string savepath = saveDir2 + filename + "_network.txt";
		std::ofstream out(savepath);
		nbn_graph.outputNBNnetwork(out);
		out.close();
		ouputNBNdata(saveDir2 + filename + "_nbn.txt", nbn_graph.get_nbn_data());
	}



}


void genenratePro(
	const std::string& proname,
	ofec::ParameterMap& params,
	std::shared_ptr<ofec::Environment>& env) {
	using namespace ofec;


	env.reset(Environment::create());
	env->recordInputParameters();
	env->initialize(0.5);
	env->setProblem(ofec::Factory<ofec::Problem>::produce(proname));
	env->problem()->inputParameters().input(params);
	env->problem()->setName(proname);
	env->problem()->recordInputParameters();
	env->initializeProblem(0.5);
}



void runAlg(const std::string& nbnSaveDir,
	const std::string& networkSaveDir) {
	using namespace ofec;
	using namespace std;
	ofec::ParameterMap params;
	std::string algName = "CMEAS-data-sampling";
	std::string proname = "RW-CEC2011-T01";
	params["problem name"] = proname;
	params["maximum evaluations"] = int(200000);
	int numRun = 160;
	//numRun = 1;

	std::shared_ptr<ofec::Environment> env;
	genenratePro(proname, params, env);
	auto pro = env->problem();
	auto rnd = make_shared<ofec::Random>(0.5);


	std::vector<std::shared_ptr<ofec::SolutionBase>> sols;
	std::vector<std::shared_ptr<ofec::SamplingData::SolutionInfo>> solInfos;
	ofec::SamplingAgorithmMutithread::runAlgMultiTask(sols, solInfos, 0.5, 0.5, proname, params, algName, params, numRun);


	std::vector<std::shared_ptr<ofec::SolutionBase>> optSols;
	std::vector<std::shared_ptr<ofec::SamplingData::SolutionInfo>> optSolInfos;
	if (pro->optimaBase()->isSolutionGiven()) {
		for (int idx(0); idx < pro->optimaBase()->numberSolutions(); ++idx) {
			optSols.emplace_back(pro->createSolution(pro->optimaBase()->solutionBase(idx)));
		}
	}



	auto filename = "Task_proName_" + proname + "_algName_" + algName;
	std::vector<ofec::SolutionBase*> solbases;

	for (auto& it : sols) {
		solbases.push_back(it.get());
	}
	for (auto& it : optSols) {
		solbases.push_back(it.get());
	}

	auto eval_fun =
		[](ofec::SolutionBase& sol, ofec::Environment* env) {
		using namespace ofec;
		sol.evaluate(env, false);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};

	UTILITY::evaluateRandomSolutionsMultiThreads(solbases, env.get(), eval_fun);
	ofec::nbn::NBNinfo nbnInfo;
	std::vector<int> mapSolId;

	CAST_CONOP(pro)->fitlerSameSolutions(solbases, nbnInfo.m_solbases, mapSolId);

	calculateNBN(nbnInfo, nbnSaveDir, networkSaveDir, filename, env.get(), rnd.get());
	{
		std::cout << "output alginfo" << std::endl;
		ofec::ParameterVariantStream paramsStream;


		paramsStream << solInfos.size();
		for (auto& it : solInfos) {
			it->toParameterVariants(paramsStream);
			//pro->solutionToParameterVariants(*it, paramsStream);
		}
		paramsStream << mapSolId;
		ofec::variants_stream::outputToFile(paramsStream, nbnSaveDir + filename + "_algInfo.txt");
	}
}


void runTask(int argc, char* argv[]) {


	// 打印参数个数和每个参数的值
	std::cout << "Number of command-line arguments: " << argc << std::endl;
	for (int i = 0; i < argc; ++i) {
		std::cout << "Argument " << i << ": " << argv[i] << std::endl;
		std::string filename = std::string(argv[i]) + "_yes";
		std::cout << "filename\t" << filename << std::endl;
	}



	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "//172.29.41.69/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "//172.29.203.176/e/DiaoYiya/code/data/ofec-data";
	//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data";
	//ofec::g_working_directory = "//172.29.203.176/e/DiaoYiya/code/data/ofec-data";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";


	auto saveDir = ofec::g_working_directory + "/paper_con_realworld/NBN_trait2/";
	auto saveDir2 = ofec::g_working_directory + "/paper_con_realworld/NBN_trait_network2/";

	std::filesystem::create_directories(saveDir);
	std::filesystem::create_directories(saveDir2);
	//std::cout << "runOneAlg\t" << argv[1] << "\t" << argv[2] << std::endl;
	//runOneAlg(saveDir, argv[1], std::stoi(argv[2]), 30);

	//std::cout << "ver5 evaluateAllAlg\t" << argv[1]  << std::endl;
	//evaluateAllAlg(saveDir, argv[1], 30);

	std::string avg = "run1";


	std::cout << "realworld algorithms and problems\t" << avg << std::endl;

	//runAlg(saveDir, saveDir2);

	//evaluateAllAlg(saveDir, "test", 3);

	//runOneAlg(saveDir, "test", 1, 30);

	//runTask();
	//runTotalTasks(std::string(argv[1]));
	//runTotalTasks("test");
}



namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run(int argc, char* argv[]) {
		registerInstance();


		runTask(argc, argv);

	}


}