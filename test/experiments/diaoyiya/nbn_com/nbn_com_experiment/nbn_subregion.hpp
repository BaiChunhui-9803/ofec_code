
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






using namespace std;





void sampleSolsAroundTSPThreadTask(
	const ofec::TravellingSalesman::SolutionType& sol,
	std::vector<std::shared_ptr<ofec::SolutionBase>>& samples,
	int neighbor_k,
	int from, int to,
	ofec::Environment* env,
	ofec::Random* rnd
) {

	auto pro = env->problem();


	int a(0), b(0);
	int dim = pro->numberVariables();

	std::shared_ptr<ofec::SolutionBase> cursol;
	for (int idx(from); idx < to; ++idx) {
		cursol.reset(pro->createSolution(sol));
		auto& cursolt = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*cursol);
		for (int idy(0); idy < neighbor_k; ++idy) {
			a = rnd->uniform.nextNonStd<int>(0, dim);
			b = rnd->uniform.nextNonStd<int>(0, dim);
			std::swap(cursolt.variable()[a], cursolt.variable()[b]);
		}

		samples[idx] = cursol;

	}
}



void sampleSolsAroundTSPMultiThread(
	const ofec::TravellingSalesman::SolutionType& sol,
	std::vector<std::shared_ptr<ofec::SolutionBase>>& samples,
	int neighbor_k,
	int numSamples,
	ofec::Environment* env,
	ofec::Random* rnd
) {

	samples.resize(numSamples);


	std::cout << "generate solutions by multithread" << std::endl;
	int num_task = std::thread::hardware_concurrency();
	int num_samples = numSamples;
	std::vector<std::thread> thrds;
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);


	std::vector<std::shared_ptr<ofec::Random>> rnds(num_task);
	for (auto& it : rnds) {
		double randomSeed(rnd->uniform.nextNonStd<double>(0.01, 1.0));
		it.reset(new ofec::Random(randomSeed));
	}

	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];

		thrds.push_back(std::thread(
			sampleSolsAroundTSPThreadTask,
			std::cref(sol), std::ref(samples), neighbor_k,
			tasks[i], tasks[i + 1], env, rnds[i].get()));
	}
	for (auto& thrd : thrds)
		thrd.join();
}







void sampleSolsRandomThreadTask(
	std::vector<std::shared_ptr<ofec::SolutionBase>>& samples,
	int from, int to,
	ofec::Environment* env,
	ofec::Random* rnd
) {
	auto pro = env->problem();
	int dim = pro->numberVariables();

	std::shared_ptr<ofec::SolutionBase> cursol;
	for (int idx(from); idx < to; ++idx) {
		cursol.reset(pro->createSolution());
		cursol->initialize(env, rnd);
		//	auto& cursolt = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*cursol);
		samples[idx] = cursol;
		cursol->evaluate(env, false);
	}
}



void sampleSolsRandomMultiThread(
	std::vector<std::shared_ptr<ofec::SolutionBase>>& samples,
	int numSamples,
	ofec::Environment* env,
	ofec::Random* rnd
) {
	samples.resize(numSamples);


	std::cout << "generate solutions by multithread" << std::endl;
	int num_task = std::thread::hardware_concurrency();
	int num_samples = numSamples;
	std::vector<std::thread> thrds;
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);


	std::vector<std::shared_ptr<ofec::Random>> rnds(num_task);
	for (auto& it : rnds) {
		double randomSeed(rnd->uniform.nextNonStd<double>(0.01, 1.0));
		it.reset(new ofec::Random(randomSeed));
	}

	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];

		thrds.push_back(std::thread(
			sampleSolsRandomThreadTask,
			std::ref(samples),
			tasks[i], tasks[i + 1], env, rnds[i].get()));
	}
	for (auto& thrd : thrds)
		thrd.join();
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


void outputToFile(ofec::ParameterVariantStream& paramsStream, const std::string& filepath) {
	std::stringstream buf;
	ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
	std::ofstream out(filepath);
	out << buf.rdbuf();
	out.close();
}

struct NBNinfo {


	std::vector<ofec::SolutionBase*> solbases;
	std::vector<int> belong;
	std::vector<double> dis2parent;
	std::vector<double> vFitness;


	void outputNBNinfo(const std::string& saveDir, const std::string& filename) {


		ofec::ParameterVariantStream paramsStream;
		paramsStream << belong;
		paramsStream << dis2parent;
		paramsStream << vFitness;
		auto filepath = saveDir + filename + "_nbnInfo.txt";
		outputToFile(paramsStream, filepath);

	}

	void outputTspSolutions(const std::string& saveDir, const std::string& filename) {
		ofec::ParameterVariantStream paramsStream;
		paramsStream << solbases.size();
		for (auto& it : solbases) {
			auto& cursol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*it);
			paramsStream << cursol.variable().vect();
		}

		outputToFile(paramsStream, saveDir + filename + "_solVariables.txt");


		//std::stringstream buf;
		//ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
		//std::cout << "buf size\t" << buf.str().size() << std::endl;
		//std::ofstream out(saveDir + filename + "_solVariables.txt");
		//out << buf.rdbuf();
		//out.close();
	}
};





void calculateNBN(NBNinfo& info, ofec::Environment* env, ofec::Random* rnd) {
	auto& solbases = info.solbases;
	nbn::HnswModel hnswModel;
	hnswModel.initialize(env, rnd, std::thread::hardware_concurrency());
	std::vector<int> solIds;
	//	new_datum->m_hnswModel2.setNumberThreads(1);
	hnswModel.addDataMutliThread(solbases, solIds);


	auto& belong = info.belong;
	auto& dis2parent = info.dis2parent;
	auto& vFitness = info.vFitness;
	vFitness.clear();

	for (auto& it : solbases) {
		vFitness.push_back(it->fitness());
	}

	nbn::HnswModelNBN model2;
	model2.copy(hnswModel);
	//model2.setNumberThreads(1);
	model2.updateFitness(vFitness);
	model2.calculateNBN(belong, dis2parent);

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
	int nbnSamples = 1e6;

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

	{
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
		NBNinfo nbnInfo;
		filterUniqueSols(solbases, nbnInfo.solbases, filterSolId, tsp_pro, rnd.get());

		UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.solbases, env.get(), eval_fun);

		for (auto& it : solIds) {
			for (auto& it2 : it) {
				for (auto& it3 : it2) {
					it3 = filterSolId[it3];
				}
			}
		}

		calculateNBN(nbnInfo, env.get(), rnd.get());
		nbnInfo.outputNBNinfo(saveDir, tspname + "_eaxRun");
		nbnInfo.outputTspSolutions(saveDir, tspname + "_eaxRun");
		{

			ofec::ParameterVariantStream paramsStream;
			paramsStream << solIds.size();
			for (const auto& it : solIds) {
				paramsStream << it.size();
				for (const auto& it2 : it) {
					paramsStream << it2;
				}
			}
			outputToFile(paramsStream, saveDir + tspname + "_eaxRun_solIds.txt");
		}

		bestSol = nbnInfo.solbases.front()->getSharedPtr();
		for (auto& it : nbnInfo.solbases) {
			if (bestSol->fitness() < it->fitness()) {
				bestSol = it->getSharedPtr();
			}
		}

	}



	{
		std::vector<std::shared_ptr<ofec::SolutionBase>> samples;
		sampleSolsRandomMultiThread(samples, numSamples, env.get(), rnd.get());
		std::vector<SolutionBase*> solbases;
		for (auto& it : samples) {
			solbases.push_back(it.get());
		}
		std::vector<int> filterSolId;
		NBNinfo nbnInfo;
		filterUniqueSols(solbases, nbnInfo.solbases, filterSolId, tsp_pro, rnd.get());
		if (nbnInfo.solbases.size() > nbnSamples)
			nbnInfo.solbases.resize(nbnSamples);

		calculateNBN(nbnInfo, env.get(), rnd.get());
		nbnInfo.outputNBNinfo(saveDir, tspname + "_randomSamples");
		nbnInfo.outputTspSolutions(saveDir, tspname + "_randomSamples");

	}

	std::vector<int> neighbors = { 100, 50,25,12,6,3 };
	auto& bestTspSol = dynamic_cast<TravellingSalesman::SolutionType&>(*bestSol);
	for (auto& neighborK : neighbors) {
		std::vector<std::shared_ptr<ofec::SolutionBase>> samples;
		sampleSolsAroundTSPMultiThread(bestTspSol, samples, neighborK, numSamples, env.get(), rnd.get());
		std::vector<SolutionBase*> solbases;
		for (auto& it : samples) {
			solbases.push_back(it.get());
		}
		std::vector<int> filterSolId;
		NBNinfo nbnInfo;
		filterUniqueSols(solbases, nbnInfo.solbases, filterSolId, tsp_pro, rnd.get());
		if (nbnInfo.solbases.size() > nbnSamples)
			nbnInfo.solbases.resize(nbnSamples);

		calculateNBN(nbnInfo, env.get(), rnd.get());
		nbnInfo.outputNBNinfo(saveDir, tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK));
		nbnInfo.outputTspSolutions(saveDir, tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK));
	}

}






void runTask() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	// 		ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";

	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_data/";

	std::filesystem::create_directories(saveDir);


	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames =
	{ "u574.tsp" , "1281.tsp" , "6702.tsp" , "2202.tsp", "6717.tsp" ,"7310.tsp", "9225.tsp", "5955.tsp", };

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




	}


}