
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


void filterUniqueSols(std::vector<std::vector<std::vector<std::shared_ptr<ofec::SolutionBase>>>>& sol,
	std::vector<std::vector<std::vector<int>>>& solIds,
	std::vector<ofec::SolutionBase*>& solbases,
	ofec::TravellingSalesman* tsp_pro, ofec::Random* rnd) {
	solIds.resize(sol.size());
	for (int idx(0); idx < solIds.size(); ++idx) {
		solIds[idx].resize(sol[idx].size());
		for (int idy(0); idy < solIds[idx].size(); ++idy) {
			solIds[idx][idy].resize(sol[idx][idy].size());
		}
	}
	//int totalSols = 0;
//	std::map<unsigned long long, size_t> solHashToId;

	ofec::TravellingSalesman::HashSolutionMap solMap;
	solMap.initialize(rnd, tsp_pro->numberVariables() + 10);




	for (int idx(0); idx < sol.size(); ++idx) {
		for (int idy(0); idy < sol[idx].size(); ++idy) {
			for (int idz(0); idz < sol[idx][idy].size(); ++idz) {
				//auto hash = tsp_pro->calHash(*sol[idx][idy][idz]);
				solIds[idx][idy][idz] = solMap.getSolId(*sol[idx][idy][idz]);
				if (solIds[idx][idy][idz] >= solbases.size()) {
					solbases.push_back(sol[idx][idy][idz].get());
				}
				//if (solHashToId.find(hash) == solHashToId.end()) {
				//	solHashToId[hash] = solbases.size();
				//	solbases.push_back(sol[idx][idy][idz].get());
				//}
				//solIds[idx][idy][idz] = solHashToId[hash];
			}
		}
	}

}



void filterUniqueSols(std::vector<std::vector<std::shared_ptr<ofec::SolutionBase>>>& sol,
	std::vector<std::vector<int>>& solIds,
	std::vector<ofec::SolutionBase*>& solbases,
	ofec::TravellingSalesman* tsp_pro, ofec::Random* rnd) {
	solIds.resize(sol.size());
	for (int idx(0); idx < solIds.size(); ++idx) {
		solIds[idx].resize(sol[idx].size());
	}
	//int totalSols = 0;
//	std::map<unsigned long long, size_t> solHashToId;

	ofec::TravellingSalesman::HashSolutionMap solMap;
	solMap.initialize(rnd, tsp_pro->numberVariables() + 10);




	for (int idx(0); idx < sol.size(); ++idx) {
		for (int idy(0); idy < sol[idx].size(); ++idy) {
			solIds[idx][idy] = solMap.getSolId(*sol[idx][idy]);
			if (solIds[idx][idy] >= solbases.size()) {
				solbases.push_back(sol[idx][idy].get());
			}
		}
	}


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


void calTaskMultiThreadInputOutputTest(const std::string& readDir,
	const std::string& saveDir,
	const std::string& tspfilename, int numRun = 30) {
	using namespace ofec;
	using namespace std;
	using namespace chrono;
	auto filename = tspfilename.substr(0, tspfilename.find_first_of("."));

	std::cout << "calculating filename\t" << filename << std::endl;

	std::shared_ptr<Environment> env;
	genenrateTSPenv(readDir, tspfilename, env);

	std::string proname = "TSP";
	auto pro = env->problem();


	//std::shared_ptr<ofec::Problem> pro(ofec::Factory<ofec::Problem>::produce(proname));
	//pro->inputParameters().input(params);
	//pro->recordInputParameters();
	////  pro->initialize(0.5);


	//auto env = std::make_shared<ofec::Environment>();
	//env->initialize();
	//env->setProblem(pro);
	//
	//env->initializeProblem(0.5);


	if (pro->numberVariables() > 2000)return;
	auto rnd = std::make_shared<ofec::Random>(0.5);
	std::string algName = "EAX-TSP-sampling";

	ofec::ParameterMap params;
	//params["problem name"] = std::string("TSP");
	params["dataFile1"] = tspfilename;
	params["dataDirectory"] = readDir;

	//params["algorithm name"] = std::string("EAX-TSP-sampling");


	auto sols = ofec::SamplingMutithread::runAlgMultiTask(0.5, 0.5, proname, params, algName, params, numRun);
	// auto sols = ofec::SamplingMutithread::runAlgMultiTask(pro.get(), algName, params, numRun);


	 //std::vector<int> popSols;
	std::vector<ofec::SolutionBase*> solbases;
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
	std::vector<std::vector<std::vector<int>>> solIds;

	filterUniqueSols(sols, solIds, solbases, tsp_pro, rnd.get());



	auto eval_fun =
		[](ofec::SolutionBase& sol, ofec::Environment* env) {
		using namespace ofec;
		sol.evaluate(env);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};
	UTILITY::evaluateRandomSolutionsMultiThreads(solbases, env.get(), eval_fun);


	std::vector<double> fitnessRuns(numRun, std::numeric_limits<double>::lowest());



	for (int idx(0); idx < fitnessRuns.size(); ++idx) {
		auto& curMax = fitnessRuns[idx];
		for (auto& itIds : solIds[idx]) {
			for (auto& itId : itIds) {
				curMax = std::max(curMax, solbases[itId]->fitness());
			}
		}
	}


	double maxFit(0);
	ofec::calMax(fitnessRuns, maxFit);

	int totalTest = 0;
	for (auto& it : fitnessRuns) {
		if (it == maxFit) {
			++totalTest;
		}
	}

	std::cout << "success Run\t" << totalTest << std::endl;


	{
		ofec::ParameterVariantStream paramsStream;
		paramsStream << solIds.size();
		for (auto& it : solIds) {
			paramsStream << it.size();
			for (auto& it2 : it) {
				paramsStream << it2;
			}
		}
		std::stringstream buf;
		ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
		std::cout << "buf size\t" << buf.str().size() << std::endl;
		std::ofstream out(saveDir + "solIds_" + tspfilename + ".txt");
		out << buf.rdbuf();
		out.close();
	}


	{
		ofec::ParameterVariantStream paramsStream;
		paramsStream << solbases.size();
		for (auto& it : solbases) {
			auto& cursol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*it);
			paramsStream << cursol.variable().vect();
		}
		std::stringstream buf;
		ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
		std::cout << "buf size\t" << buf.str().size() << std::endl;
		std::ofstream out(saveDir + "solVariables_" + tspfilename + ".txt");
		out << buf.rdbuf();
		out.close();
	}




	{
		ofec::ParameterVariantStream paramsStream;
		paramsStream << fitnessRuns;
		std::stringstream buf;
		ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
		std::cout << "buf size\t" << buf.str().size() << std::endl;
		std::ofstream out(saveDir + "algorithmBestValue_" + tspfilename + ".txt");
		out << buf.rdbuf();
		out.close();
	}



}


struct LKH_RunInfo {
	std::vector<unsigned> m_seeds;
	std::vector<std::vector<std::vector<int>>> m_sols;

	int m_curId = 0;
	std::mutex m_mutex;


	void resetSize(int numRun) {
		m_seeds.resize(numRun);
		m_sols.resize(numRun);
	}

};


void runLKH(const std::string& readDir,
	const std::string& tspfilename,
	LKH_RunInfo& curInfo) {

	int curRunId = 0;
	while (true) {
		{
			std::unique_lock lock(curInfo.m_mutex);
			if (curInfo.m_curId == curInfo.m_seeds.size())return;
			curRunId = curInfo.m_curId++;
		}


		auto& Seed = curInfo.m_seeds[curRunId];
		auto& totalSols = curInfo.m_sols[curRunId];
		LKH::LKHAlg lkh_alg;

		lkh_alg.SRandom(Seed);
		lkh_alg.run(readDir, tspfilename, totalSols);

		for (auto& it : totalSols) {
			for (auto& it2 : it) {
				--it2;
			}
		}
	}




}


void runLKHmultiThread(
	const std::string& readDir,
	const std::string& saveDir,
	const std::string& tspfilename) {
	using namespace ofec;
	int numRuns(30);
	//numRuns = 1;
	LKH_RunInfo totalInfo;
	totalInfo.resetSize(numRuns);

	unsigned initSeed = 1;
	unsigned SeedStep = 1e4;

	for (auto& it : totalInfo.m_seeds) {
		it = initSeed;
		initSeed += SeedStep;
	}

	int num_task = std::thread::hardware_concurrency();
	std::vector<std::thread> thrds;
	for (size_t i = 0; i < num_task; ++i) {
		thrds.push_back(std::thread(
			runLKH,
			std::cref(readDir), std::cref(tspfilename), std::ref(totalInfo)));
	}
	for (auto& thrd : thrds)
		thrd.join();


	std::vector<std::vector<std::shared_ptr<ofec::SolutionBase>>> totalSols;

	std::shared_ptr<Environment> env;
	genenrateTSPenv(readDir, tspfilename, env);
	auto pro = env->problem();

	auto rnd = std::make_shared<Random>(0.5);

	totalSols.resize(numRuns);
	for (int idx(0); idx < totalSols.size(); ++idx) {
		auto& curAlg = totalSols[idx];
		auto& curAlgSol = totalInfo.m_sols[idx];
		curAlg.resize(curAlgSol.size());
		for (int idy(0); idy < curAlg.size(); ++idy) {
			auto& curIter = curAlg[idy];
			auto& curIterSol = curAlgSol[idy];

			curIter.reset(pro->createSolution());
			auto& tspSol = dynamic_cast<TravellingSalesman::SolutionType&>(*curIter);
			tspSol.variable().vect() = curIterSol;
		}
	}


	//std::vector<int> popSols;
	std::vector<ofec::SolutionBase*> solbases;
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
	std::vector<std::vector<int>> solIds;

	filterUniqueSols(totalSols, solIds, solbases, tsp_pro, rnd.get());



	auto eval_fun =
		[](ofec::SolutionBase& sol, ofec::Environment* env) {
		using namespace ofec;
		sol.evaluate(env);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};
	UTILITY::evaluateRandomSolutionsMultiThreads(solbases, env.get(), eval_fun);


	std::vector<double> fitnessRuns(numRuns, std::numeric_limits<double>::lowest());

	for (int idx(0); idx < fitnessRuns.size(); ++idx) {
		auto& curMax = fitnessRuns[idx];
		for (auto& itId : solIds[idx]) {
			curMax = std::max(curMax, solbases[itId]->fitness());

		}
	}

	double maxFit(0);
	ofec::calMax(fitnessRuns, maxFit);

	int totalTest = 0;
	for (auto& it : fitnessRuns) {
		if (it == maxFit) {
			++totalTest;
		}
	}

	std::cout << "success Run\t" << totalTest << std::endl;

	{
		ofec::ParameterVariantStream paramsStream;
		paramsStream << solIds.size();
		for (auto& it : solIds) {
			paramsStream << it;
		}
		std::stringstream buf;
		ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
		std::cout << "buf size\t" << buf.str().size() << std::endl;
		std::ofstream out(saveDir + "solIds_" + tspfilename + ".txt");
		out << buf.rdbuf();
		out.close();
	}


	{
		ofec::ParameterVariantStream paramsStream;
		paramsStream << solbases.size();
		for (auto& it : solbases) {
			auto& cursol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*it);
			paramsStream << cursol.variable().vect();
		}
		std::stringstream buf;
		ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
		std::cout << "buf size\t" << buf.str().size() << std::endl;
		std::ofstream out(saveDir + "solVariables_" + tspfilename + ".txt");
		out << buf.rdbuf();
		out.close();
	}



	{
		ofec::ParameterVariantStream paramsStream;
		paramsStream << fitnessRuns;
		std::stringstream buf;
		ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
		std::cout << "buf size\t" << buf.str().size() << std::endl;
		std::ofstream out(saveDir + "algorithmBestValue_" + tspfilename + ".txt");
		out << buf.rdbuf();
		out.close();
	}




}





void readDir(const std::string& readDir, std::vector<std::string>& filenames) {

	for (const auto& entry : std::filesystem::directory_iterator(readDir)) {
		std::cout << entry.path().filename() << std::endl;
	}

}



namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run() {

		using namespace ofec;
		using namespace std;

		registerInstance();
		ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
		//	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";


		std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
		saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
		saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_data/";

		std::filesystem::create_directories(saveDir);


		std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
		std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";



		//runLKHmultiThread(dir1, saveDir, "2202.tsp");
		//runLKHmultiThread(dir2, saveDir, "u574.tsp");

		std::vector<std::string> filenames;

		readDir(dir1, filenames);

		//calTaskMultiThreadInputOutputTest(dir1, saveDir, "2202.tsp");
		//calTaskMultiThreadInputOutputTest(dir2, saveDir, "u574.tsp");




	//	runLKH(dir2, saveDir, "GR666.tsp", 1);


		//runTask(saveDir, dir2, "GR666.tsp");
		//runTask(saveDir, dir1, "2202.tsp");


	//	runLonGen(dir2, saveDir, "GR666.tsp");
	//	runLonGen(dir1, saveDir, "2202.tsp");


		//inputLonInfo(saveDir, "GR666.tsp");

	}


}