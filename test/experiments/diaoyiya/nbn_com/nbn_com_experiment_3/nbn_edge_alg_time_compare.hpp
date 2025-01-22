
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
#include "../utility/nbn_visualization/nbn_calculator/nbn_modifiy_equal.h"


//#include "../utility/nbn_visualization/calculate_algorithm/"
#include "../utility/nbn_visualization/calculate_algorithm/nbn_edge_multithread_division.h"
#include "../core/parameter/variants_to_stream.h"

#include "../utility/nbn_visualization/calculate_algorithm/nbn_iteration_multithread.h"



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

struct NBNAlgInfo {
	int m_iteration = 0;
	int m_numLoop = 0;
	double m_maxChange = 0;
	double m_dev = 0;
	double m_millisecond = 0;

	static void outputHead(std::ostream& out) {
		out << "iteration\tloop\tmax\tdev\tseconds" << std::endl;
	}
	void output(std::ostream& out) {
		out << m_iteration << "\t" << m_numLoop << "\t" << m_maxChange << "\t" << m_dev << "\t" << m_millisecond << std::endl;
	}

	void updateInfo(const std::vector<double>& dis1, const std::vector<double>& dis2) {
		m_dev = 0;
		m_maxChange = 0;
		for (int idx(0); idx < dis1.size(); ++idx) {
			m_maxChange = std::max(m_maxChange, std::abs(dis1[idx] - dis2[idx]));
			m_dev += (dis1[idx] - dis2[idx]) * (dis1[idx] - dis2[idx]);
		}

		m_dev = sqrt(m_dev / dis1.size());

	}
};


void solDataToNBN(
	const std::string& outputDir,
	const std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
	const std::string& filename,
	ofec::Environment* env,
	ofec::Random* rnd
	//,
	//std::vector<NBNAlgInfo> & algInfo
) {
	using namespace ofec;


	std::shared_ptr<ofec::NBN_EdgeMultiThreadDivision> m_division(new ofec::NBN_EdgeMultiThreadDivision);


	auto eval_fun =
		[](ofec::SolutionBase& sol, Environment* env) {
		using namespace ofec;
		sol.evaluate(env, false);
		auto pro = env->problem();
		ofec::Real pos = (pro->optimizeMode(0) == ofec::OptimizeMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};
	{
		auto& division(m_division);
		division->setMaxSampleSize(0);
		division->setMultiThread(true);
		int numThread = std::thread::hardware_concurrency();
		division->setNumberThread(numThread);
		std::cout << "number of thread" << numThread << std::endl;
		division->initialize(env->getSharedPtr(), rnd->getSharedPtr(), eval_fun, true);
		auto& edge_division = dynamic_cast<ofec::NBN_EdgeMultiThreadDivision&>(*division);
		int numLoop = 1e2 + 5;
		//edge_division.setNumberLoop(std::ceil(double(numLoop) / double(numThread)));
		edge_division.setNumberLoop(1);
		edge_division.updateSols(sols);
		//int from_size = sols.size();
		//edge_division.generateSols(sols.size() + 1e3);
		//int to_size = edge_division.size();
		//edge_division.filterSameSolutions(from_size, to_size);

		auto start = std::chrono::system_clock::now();
		edge_division.updateNetwork();
		auto end = std::chrono::system_clock::now();
		//std::chrono::duration<float> difference = end - start;
	//	double milliseconds = difference.count();

		//auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(difference).count();
		std::chrono::duration<double> diffSeconds = end - start; // 单位为秒
		auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(diffSeconds).count();
		std::vector<double> beforeDis2parent(sols.size(), 1);
		std::vector<int> belong;
		std::vector<double> dis2parent;
		edge_division.getResult(belong, dis2parent);

		std::ofstream algInfoOut(outputDir + filename + "_algInfo.txt");
		NBNAlgInfo algInfo;
		NBNAlgInfo::outputHead(algInfoOut);
		algInfo.m_iteration = 0;
		algInfo.m_numLoop = numThread;
		algInfo.m_millisecond = differenceInSeconds;
		algInfo.updateInfo(beforeDis2parent, dis2parent);
		algInfo.output(algInfoOut);
		algInfo.output(std::cout);





		while (true) {
			swap(dis2parent, beforeDis2parent);
			edge_division.updateDivisionSubRegionMultiThread();
			//end = std::chrono::system_clock::now();
			auto end = std::chrono::system_clock::now();
			//std::chrono::duration<float> difference = end - start;
		//	double milliseconds = difference.count();

			//auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(difference).count();
			std::chrono::duration<double> diffSeconds = end - start; // 单位为秒
			auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(diffSeconds).count();
			edge_division.getResult(belong, dis2parent);

			algInfo.m_iteration++;
			algInfo.m_millisecond = differenceInSeconds;
			algInfo.updateInfo(beforeDis2parent, dis2parent);
			algInfo.output(algInfoOut);
			algInfo.output(std::cout);

			if (algInfo.m_maxChange < 0.01 && algInfo.m_dev < 0.1) break;
		}

		{
			std::vector<int> belongIds(belong.size());
			std::vector<double> fitness(sols.size());
			for (int idx(0); idx < belongIds.size(); ++idx) {
				belongIds[idx] = idx;
				fitness[idx] = sols[idx]->fitness();
			}


			std::ofstream nbnOut(outputDir + filename + "_nbn.txt");
			ofec::outputNBNdata(nbnOut, belongIds, belong, dis2parent, fitness);
			nbnOut.close();
		}
		//edge_division.updateBestSolSub(1e3);

		////	swap(dis2parent, beforeDis2parent);
		//std::ofstream nbnOut(filepath + ".nbn");
		//std::vector<double> fitness(sols.size());
		//for (int idx(0); idx < fitness.size(); ++idx) {
		//	fitness[idx] = sols[idx]->fitness();
		//}
		//outputNBNdata(nbnOut, solId, belong, dis2parent, fitness);
		//nbnOut.close();


	}
}




void solDataToNBN_subRegion(
	const std::string& outputDir,
	const std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
	std::shared_ptr<ofec::SolutionBase>& centerSol,
	const std::string& filename,
	ofec::Environment* env,
	ofec::Random* rnd
	//,
	//std::vector<NBNAlgInfo> & algInfo
) {
	using namespace ofec;

	//	algInfo.clear();
	std::cout << "running\t" << filename << std::endl;



	//ofec::ParamMap params;
	//params["problem name"] = std::string("TSP");
	//params["dataFile1"] = tspname;
	////	params["division type"] = static_cast<int>(ofec::NBN_VisualizationData::NBN_Division_Type::kEdgeDivision);
	//	//  params["continous sample size"] = int(2e3);
	////params["continous sample size"] = int(2e6 + 2e3);
	////params["init sample size"] = int(2e6);
	//params["add optimal"] = false;


	//auto param= ADD_PARAM(params);
	//auto pro = ADD_PRO(param, 0.1);
	//pro->initialize();

	//auto rnd = std::make_shared<Random>(0.5);



	std::shared_ptr<ofec::NBN_EdgeMultiThreadDivision> m_division(new ofec::NBN_EdgeMultiThreadDivision);





	auto eval_fun =
		[](ofec::SolutionBase& sol, Environment* env) {
		using namespace ofec;
		sol.evaluate(env, false);
		auto pro = env->problem();
		ofec::Real pos = (pro->optimizeMode(0) == ofec::OptimizeMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};
	{
		auto& division(m_division);

		division->setMaxSampleSize(0);
		division->setMultiThread(true);
		int numThread = std::thread::hardware_concurrency();
		division->setNumberThread(numThread);
		std::cout << "number of thread" << numThread << std::endl;
		division->initialize(env->getSharedPtr(), rnd->getSharedPtr(), eval_fun, true);
		auto& edge_division = dynamic_cast<ofec::NBN_EdgeMultiThreadDivision&>(*division);
		int numLoop = 1e2 + 5;
		//edge_division.setNumberLoop(std::ceil(double(numLoop) / double(numThread)));
		edge_division.setNumberLoop(1);
		edge_division.updateSols(sols);
		//int from_size = sols.size();
		//edge_division.generateSols(sols.size() + 1e3);
		//int to_size = edge_division.size();
		//edge_division.filterSameSolutions(from_size, to_size);

		auto start = std::chrono::system_clock::now();
		//edge_division.updateNetwork();
		auto end = std::chrono::system_clock::now();
		//std::chrono::duration<float> difference = end - start;
	//	double milliseconds = difference.count();

		//auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(difference).count();
		std::chrono::duration<double> diffSeconds = end - start; // 单位为秒
		auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(diffSeconds).count();

		std::vector<double> beforeDis2parent(sols.size(), 1);
		std::vector<int> belong;
		std::vector<double> dis2parent;
		edge_division.getResult(belong, dis2parent);



		std::ofstream algInfoOut(outputDir + filename + "_algInfo.txt");
		NBNAlgInfo algInfo;
		NBNAlgInfo::outputHead(algInfoOut);

		algInfo.m_iteration = 0;
		algInfo.m_numLoop = numThread;
		algInfo.m_millisecond = differenceInSeconds;
		algInfo.updateInfo(beforeDis2parent, dis2parent);
		swap(dis2parent, beforeDis2parent);
		algInfo.output(algInfoOut);
		algInfo.output(std::cout);

		edge_division.setCenterSol(centerSol);
		edge_division.updateBestSolSub(1e3);
		//	edge_division.updateSolsDis2Center();

		int minLoop = 5;
		while (true) {
			minLoop--;
			edge_division.updateDivisionSubRegionBandEdgeMultiThread();
			auto end = std::chrono::system_clock::now();
			//std::chrono::duration<float> difference = end - start;
		//	double milliseconds = difference.count();

			//auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(difference).count();
			std::chrono::duration<double> diffSeconds = end - start; // 单位为秒
			auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(diffSeconds).count();
			edge_division.getResult(belong, dis2parent);

			algInfo.m_iteration++;
			algInfo.m_millisecond = differenceInSeconds;
			algInfo.updateInfo(beforeDis2parent, dis2parent);
			swap(dis2parent, beforeDis2parent);
			algInfo.output(algInfoOut);
			algInfo.output(std::cout);

			if (minLoop < 0 && algInfo.m_maxChange < 0.01 && algInfo.m_dev < 0.1) break;
		}

		{
			std::vector<int> belongIds(belong.size());
			std::vector<double> fitness(sols.size());
			for (int idx(0); idx < belongIds.size(); ++idx) {
				belongIds[idx] = idx;
				fitness[idx] = sols[idx]->fitness();
			}


			std::ofstream nbnOut(outputDir + filename + "_nbn.txt");
			ofec::outputNBNdata(nbnOut, belongIds, belong, dis2parent, fitness);
			nbnOut.close();
		}

		//dis2parent = beforeDis2parent;

		//std::ofstream nbnOut(filepath + ".nbn");
		//std::vector<double> fitness(sols.size());
		//for (int idx(0); idx < fitness.size(); ++idx) {
		//	fitness[idx] = sols[idx]->fitness();
		//}
		//outputNBNdata(nbnOut, solId, belong, dis2parent, fitness);
		//nbnOut.close();


	}
}

void solDataToNBN_itererationMultithread(
	const std::string& outputDir,
	const std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
	std::shared_ptr<ofec::SolutionBase>& centerSol,
	const std::string& filename,
	ofec::Environment* env,
	ofec::Random* rnd
	//,
	//std::vector<NBNAlgInfo> & algInfo
) {
	std::vector<ofec::SolutionBase*> solbases;
	for (auto& it : sols) {
		solbases.push_back(it.get());

	}
	std::vector<int> belong;
	std::vector<double> dis2parent;


	std::ofstream algInfoOut(outputDir + filename + "_algInfo.txt");
	NBNAlgInfo algInfo;
	NBNAlgInfo::outputHead(algInfoOut);
	algInfo.output(std::cout);



	auto start = std::chrono::system_clock::now();


	ofec::nbn::calculateNBNiterationMultithread(
		solbases, belong, dis2parent,
		std::thread::hardware_concurrency(),
		env, rnd);

	//edge_division.updateNetwork();
	auto end = std::chrono::system_clock::now();
	//std::chrono::duration<float> difference = end - start;
//	double milliseconds = difference.count();

	//auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(difference).count();
	std::chrono::duration<double> diffSeconds = end - start; // 单位为秒
	auto differenceInSeconds = std::chrono::duration_cast<std::chrono::seconds>(diffSeconds).count();


	algInfo.m_iteration = 0;
	algInfo.m_numLoop = std::thread::hardware_concurrency();
	algInfo.m_millisecond = differenceInSeconds;
	algInfo.output(algInfoOut);
	algInfo.output(std::cout);
}


void calTaskOneTask(
	const std::string& saveDir,
	const std::string& tspfilename,
	ofec::SolutionBase* bestSol,
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
	auto sampleRand = std::make_shared < ofec::Random>(0.5);
	if (neighborK < 0) {
		filename = tspname + "_randomSamples";

		std::string filepath = { nbn::NBNinfo::getSolutionsFileName(saveDir, filename) };
		bool filepathExists = std::filesystem::exists(filepath);
		if (!filepathExists) {

			std::cout << "generate solutions by multithread " << filename << std::endl;
			ofec::nbn::tsp::sampleSolsRandomMultiThread(samples, numSamples, env, sampleRand.get());
		}
	}
	else {
		filename = tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK);

		{
			auto& bestTspSol = dynamic_cast<TravellingSalesman::SolutionType&>(*bestSol);
			std::cout << "generate solutions by multithread sampling neighborK\t" << neighborK << std::endl;


			std::string filepath = { NBNinfo::getSolutionsFileName(saveDir,tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK)) };
			bool filepathExists = std::filesystem::exists(filepath);

			if (!filepathExists) {
				//		std::vector<std::shared_ptr<ofec::SolutionBase>> samples;
				nbn::tsp::sampleSolsAroundTSPMultiThread(bestTspSol, samples, neighborK, numSamples, env, sampleRand.get());
			}
		}


	}


	if (!samples.empty()) {
		std::vector<SolutionBase*> solbases;
		for (auto& it : samples) {
			solbases.push_back(it.get());
		}

		NBNinfo nbnInfo;
		std::cout << "filterUniqueSols solutions by multithread" << std::endl;
		{
			std::vector<int> filterSolId;
			nbn::tsp::filterUniqueSols(solbases, nbnInfo.m_solbases, filterSolId, tsp_pro, sampleRand.get());
		}
		nbn::tsp::sortEachSolutionX(nbnInfo.m_solbases);
		UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.m_solbases, env, eval_fun);


		if (nbnInfo.m_solbases.size() > nbnSamples)
			nbnInfo.m_solbases.resize(nbnSamples);
		std::vector<std::shared_ptr<ofec::SolutionBase>> sols;
		for (auto& it : nbnInfo.m_solbases) {
			sols.push_back(it->getSharedPtr());
		}


		auto centerSol = bestSol->getSharedPtr();
		if (neighborK < 0) {
			//solDataToNBN(saveDir, sols, filename + "_random", env, rnd);
			solDataToNBN_itererationMultithread(saveDir, sols, centerSol, filename + "_iteration", env, rnd);
		}
		else {
			solDataToNBN_itererationMultithread(saveDir, sols, centerSol, filename + "_iteration", env, rnd);
		}
	}

}




void calTaskTotal(const std::string& readDir,
	const std::string& saveDir,
	const std::string& tspfilename,
	const std::vector<int>& tasks,
	int nbnSamples = 1e6
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



	//{
	//	bestSol .reset(env->problem()->createSolution());
	//	//std::ofstream out(saveDir + "bestSolutions.txt");
	//	ofec::ParameterVariantStream stream;
	//	ofec::variants_stream::inputFromFile(stream, saveDir + tspfilename + "_bestSolutions.txt");
	//	env->problem()->parameterVariantsToSolution(stream ,*bestSol);

	//	//out << stream;
	////	out.close();
	//}

	std::vector<double> div_samples = { 1e3,1e4,1e5,1e6 };

	{
		for (auto& curSample : div_samples) {
			auto filename = tspname + "numSamples_" + std::to_string(int(curSample));
			for (auto& curTask : tasks) {
				calTaskOneTask(
					saveDir, filename, bestSol.get(), curTask, curSample, env.get(), rnd.get());
			}
		}
	}


}


void runTask2() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
//	ofec::g_working_directory= "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";
//	ofec::g_working_directory = "//172.24.24.151/e/DiaoYiya/code/data/ofec-data";

	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp2/nbn_data_eax_lkh_compare/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_tsp_algorithmTime/";
	std::filesystem::create_directories(saveDir);

	//auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread_network/";
	//std::filesystem::create_directories(saveDir2);

	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames /*=
	{ "2202.tsp", "u574.tsp" ,"5955.tsp", "1281.tsp" ,"6702.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", }*/;

	filenames = { "u574.tsp",  "5955.tsp" ,   "2202.tsp", };

	std::vector<std::pair<std::string, std::vector<int>>> tasks;


	//	tasks.push_back({ "u574.tsp" , std::vector<int>({ -1,,12,6,3*/}) });
	//tasks.push_back({ "u574.tsp" , std::vector<int>({ -1,12,3}) });
	//{3,6,12,}

	tasks.push_back({ "5955.tsp" , std::vector<int>({ 100, 50,25,12,3,-1 }) });
	tasks.push_back({ "u574.tsp" , std::vector<int>({ 100, 50,25,12,3,-1 }) });
	tasks.push_back({ "2202.tsp" , std::vector<int>({ 100, 50,25,3,-1 }) });


	std::reverse(tasks.begin(), tasks.end());
	std::vector<int> neighbors = { -1, 100, 50,25,12,6,3 };
	//	neighbors = { -1 };
	for (auto& it : tasks) {
		it.second = neighbors;
		std::reverse(it.second.begin(), it.second.end());

	}

	std::cout << "calculate nbn alg time " << std::endl;

	for (auto& it : tasks) {
		if (it.first == "u574.tsp") {
			calTaskTotal(dir2, saveDir, it.first, it.second);
			//calTask(dir2, saveDir, saveDir2, it);
		}
		else {
			calTaskTotal(dir1, saveDir, it.first, it.second);
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