
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
#include "../utility/matlab/matlab_utility.h"

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
#include "../utility/nbn_visualization/calculate_algorithm/nbn_iteration_multithread.h"
#include "../utility/nbn_visualization/nbn_fla/tsp_related/nbn_fla_tsp.h"
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


	sortedIds.resize(solSize);


	auto oldDis = m_nbnInfo.m_dis2parent;

	ofec::nbn::updateNearsetBetterSolutions(
		sortedIds, m_nbnInfo.m_solbases, m_nbnInfo.m_belong, m_nbnInfo.m_dis2parent, m_nbnInfo.m_vFitness, env, rnd);
	//	m_nbnInfo, updateIds, from, to, sortedIds, env, rnd);env, rnd


	/*
	const std::vector<int>& solIds,
			const std::vector<ofec::SolutionBase*>& sols,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			const std::vector<double>& fitness,
			ofec::Environment* env,
			ofec::Random* rnd*/
	double improveAccuracy(0);

	double maxImprove(0);

	for (int idx(0); idx < m_nbnInfo.m_dis2parent.size(); ++idx) {
		if (m_nbnInfo.m_dis2parent[idx] < oldDis[idx]) {
			++improveAccuracy;
			maxImprove = std::max(maxImprove, oldDis[idx] - m_nbnInfo.m_dis2parent[idx]);
		}

		//	if (updateIds[idx])
	}
	std::cout << "improveAccuracy\t" << improveAccuracy << "\tMax Imrprove Dis\t" << maxImprove << std::endl;


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

	//udpateNBNinfoIteration(nbnInfo,
	//	1e3, env, rnd);


	nbnInfo.outputNBNinfo(saveDir, filename);
	nbnInfo.outputVecSolutions<int>(saveDir, filename, env);


	{

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



void readNBN(
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
	//	nbnInfo.calculateNBNhnswEqualBetterRandom(env, rnd);

		//udpateNBNinfoIteration(nbnInfo,
		//	1e3, env, rnd);


	nbnInfo.inputNBNinfo(saveDir, filename);
	//nbnInfo.inputVecSolutions<int>(saveDir, filename, env);



}




void calculateBelongOneTask(std::vector<ofec::SolutionBase*>& centerSols,
	ofec::nbn::NBNinfo& mergeNBNinfo,
	std::vector<int>& funnelIds,
	std::vector<double>& funnelToDis,
	int from, int to, ofec::Environment* env) {
	for (int idx(from); idx < to; ++idx) {
		auto cursol = mergeNBNinfo.m_solbases[idx];
		double maxDis = std::numeric_limits<double>::max();
		for (int idOpt(0); idOpt < centerSols.size(); ++idOpt) {
			double curdis = cursol->variableDistance(*centerSols[idOpt], env);
			if (curdis == maxDis) {
				funnelIds[idx] = -1;
			}
			else if (curdis < maxDis) {
				funnelIds[idx] = idOpt;
				maxDis = curdis;
			}
		}
	}
}

void calculateBelong(std::vector<ofec::SolutionBase*>& centerSols,
	ofec::nbn::NBNinfo& mergeNBNinfo, std::vector<int>& funnelIds,
	std::vector<double>& funnelToDis,
	ofec::Environment* env
) {
	funnelIds.resize(mergeNBNinfo.m_solbases.size());
	std::fill(funnelIds.begin(), funnelIds.end(), -1);

	int num_task = std::thread::hardware_concurrency();
	int num_samples = mergeNBNinfo.m_solbases.size();
	std::vector<std::thread> thrds;
	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);
	std::pair<int, int> from_to;


	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];
		thrds.push_back(std::thread(
			calculateBelongOneTask,
			std::ref(centerSols), std::ref(mergeNBNinfo), std::ref(funnelIds), std::ref(funnelToDis), from_to.first, from_to.second, env));
	}

	for (auto& thrd : thrds)
		thrd.join();

	//	calculateBelongOneTask(centerSols, mergeNBNinfo, funnelIds, 0,funnelIds.size(), env);

}



void calTask(const std::string& readDir,
	const std::string& traitDir,
	const std::string& saveDir,
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

	ofec::nbn::NBNinfo optNBNinfo;
	{
		auto filename = tspname + "_filterOptSols";
		optNBNinfo.inputVecSolutions<int>(traitDir, filename, env.get());
		optNBNinfo.inputNBNinfo(traitDir, filename);
	}



	std::vector<std::string> dataFilenames = { saveDir + "optima_data/", saveDir + "secondFurther_data/", saveDir + "merge_data/" };
	std::vector<std::string> networkFilenames = { saveDir + "optima_network/", saveDir + "secondFurther_network/", saveDir + "merge_network/" };

	for (auto& it : dataFilenames) {
		std::filesystem::create_directories(it);
	}

	for (auto& it : networkFilenames) {
		std::filesystem::create_directories(it);
	}


	int bestSolId = 0;
	{
		std::vector<int> solIds;
		ofec::nbn::getBestSolIds(optNBNinfo, solIds);
		bestSolId = solIds.front();
		//centerSols.push_back(optNBNinfo.m_solbases[solIds.front()]);
	}


	std::vector<int> neiKs = { 8 };
	std::vector<double> numSamplesSize = { 1e6 };
	//outputNBNdata(std::cout, optNBNinfo.m_solIds, optNBNinfo.m_belong, optNBNinfo.m_dis2parent, optNBNinfo.m_vFitness);;

	std::pair<double, double> fit_range = { std::numeric_limits<double>::max(),std::numeric_limits<double>::lowest() };
	for (int filterSolId(0); filterSolId < optNBNinfo.m_belong.size(); ++filterSolId) {

		if (filterSolId == bestSolId)continue;
		std::vector<SolutionBase*> centerSols;
		{
			//std::vector<int> solIds;
			//ofec::nbn::getBestSolIds(optNBNinfo, solIds);
			centerSols.push_back(optNBNinfo.m_solbases[bestSolId]);
		}
		{
			//std::vector<int> solIds;
			//ofec::nbn::getSecondFurtherSolIds(optNBNinfo, solIds);
			centerSols.push_back(optNBNinfo.m_solbases[filterSolId]);
		}

		{
			std::vector<int> subsetIds = { bestSolId, filterSolId };
			ofec::nbn::NBNinfo centerNBN;
			centerNBN.subset(optNBNinfo, subsetIds);
			//	outputNBNdata(std::cout, centerNBN.m_solIds, centerNBN.m_belong, centerNBN.m_dis2parent, centerNBN.m_vFitness);;


		}


		std::vector<NBNinfo> nbnInfos(centerSols.size());
		//int neighborK = 40;


		for (auto& curSample : numSamplesSize) {
			numSamples = curSample;

			for (auto& neighborK : neiKs) {


				auto filename = tspname + "_randomSamples" + "_solId_" + std::to_string(filterSolId) + "_neighborK_" + std::to_string(neighborK) + "_numSamples_" + std::to_string(numSamples);


				NBNinfo mergeNBNinfo;

				{
					readNBN(mergeNBNinfo, dataFilenames.back(), networkFilenames.back(), filename, env.get(), rnd.get());
				}
				{
					double maxValue(0), minValue(0);
					calMax(mergeNBNinfo.m_vFitness, maxValue);
					calMin(mergeNBNinfo.m_vFitness, minValue);

					fit_range.first = std::min(minValue, fit_range.first);
					fit_range.second = std::max(maxValue, fit_range.second);

				}
			}
		}
	}


	for (int filterSolId(0); filterSolId < optNBNinfo.m_belong.size(); ++filterSolId) {

		if (filterSolId == bestSolId)continue;
		std::vector<SolutionBase*> centerSols;
		{
			//std::vector<int> solIds;
			//ofec::nbn::getBestSolIds(optNBNinfo, solIds);
			centerSols.push_back(optNBNinfo.m_solbases[bestSolId]);
		}
		{
			//std::vector<int> solIds;
			//ofec::nbn::getSecondFurtherSolIds(optNBNinfo, solIds);
			centerSols.push_back(optNBNinfo.m_solbases[filterSolId]);
		}

		{
			std::vector<int> subsetIds = { bestSolId, filterSolId };
			ofec::nbn::NBNinfo centerNBN;
			centerNBN.subset(optNBNinfo, subsetIds);
			//	outputNBNdata(std::cout, centerNBN.m_solIds, centerNBN.m_belong, centerNBN.m_dis2parent, centerNBN.m_vFitness);;


		}


		std::vector<NBNinfo> nbnInfos(centerSols.size());
		//int neighborK = 40;
		//std::vector<int> neiKs = { 8 };


		//std::vector<double> numSamplesSize = { 1e6 };

		for (auto& curSample : numSamplesSize) {
			numSamples = curSample;

			for (auto& neighborK : neiKs) {


				auto filename = tspname + "_randomSamples" + "_solId_" + std::to_string(filterSolId) + "_neighborK_" + std::to_string(neighborK) + "_numSamples_" + std::to_string(numSamples);

				std::cout << filename << std::endl;
				//for (int id(0); id < centerSols.size(); ++id) {

				//	std::cout << dataFilenames[id] << std::endl;

				//	auto& bestTspSol = dynamic_cast<TravellingSalesman::SolutionType&>(*centerSols[id]);
				//	std::vector<std::shared_ptr<ofec::SolutionBase>> samples;
				//	nbn::tsp::sampleSolsAroundTSPMultiThread_kSwap(bestTspSol, samples, neighborK, numSamples, env.get(), rnd.get());
				//	std::vector<SolutionBase*> solbases;
				//	for (auto& it : samples) {
				//		solbases.push_back(it.get());
				//	}

				//	auto& nbnInfo = nbnInfos[id];
				//	{
				//		std::vector<int> filterSolId;
				//		filterUniqueSols(solbases, nbnInfo.m_solbases, filterSolId, tsp_pro, rnd.get());
				//		//nbn::tsp::sortEachSolutionX(nbnInfo.m_solbases);
				//	}

				//	if (nbnInfo.m_solbases.size() > nbnSamples)
				//		nbnInfo.m_solbases.resize(nbnSamples);


				//	{
				//		nbnInfo.m_sols.clear();
				//		for (auto& it : nbnInfo.m_solbases) {
				//			nbnInfo.m_sols.push_back(it->getSharedPtr());
				//		}
				//		nbnInfo.m_solbases.clear();
				//		for (auto& it : nbnInfo.m_sols) {
				//			nbnInfo.m_solbases.push_back(it.get());
				//		}
				//	}

				//	//nbn::tsp::sortEachSolutionX(nbnInfo.m_solbases);




				//	UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.m_solbases, env.get(), eval_fun);

				//	calculateNBN(nbnInfo, dataFilenames[id], networkFilenames[id], filename, env.get(), rnd.get());
				//}

				NBNinfo mergeNBNinfo;

				{
					readNBN(mergeNBNinfo, dataFilenames.back(), networkFilenames.back(), filename, env.get(), rnd.get());
				}

				{



					std::vector<int> funnelIds;

					auto filedir = saveDir + "funnelInfo/";
					std::filesystem::create_directories(filedir);

					{
						ofec::ParameterVariantStream paramsStream;
						variants_stream::inputFromFile(paramsStream, filedir + filename + "_funnelInfo.txt");
						paramsStream >> funnelIds;
						//variants_stream::outputToFile(paramsStream, filedir + filename + "_funnelInfo.txt");
					}


					for (int nodeId(0); nodeId < mergeNBNinfo.m_solIds.size(); ++nodeId) {
						if (mergeNBNinfo.m_belong[nodeId] != -1 && mergeNBNinfo.m_belong[nodeId] != nodeId) {
							if (funnelIds[mergeNBNinfo.m_belong[nodeId]] != funnelIds[nodeId]) {
								std::cout << "tspname\t" << tspname << "\tSolId\t" << filterSolId << "\tBridgeDis\t" << mergeNBNinfo.m_dis2parent[nodeId] << std::endl;
							}
						}
					}



					std::vector<std::vector<int>> belongIds(2);
					std::vector<std::vector<double>> belongIdFit(2);
					auto norFit = mergeNBNinfo.m_vFitness;
					//ofec::dataNormalize(norFit);
					ofec::dataNormalizeInBound(norFit, fit_range);


					for (int idx(0); idx < funnelIds.size(); ++idx) {
						if (funnelIds[idx] >= 0) {
							belongIds[funnelIds[idx]].push_back(idx);
							belongIdFit[funnelIds[idx]].push_back(norFit[idx]);
						}
					}


					for (int idx(0); idx < belongIdFit.size(); ++idx) {
						double mean(0), stdValue(0);
						calMeanAndStd(belongIdFit[idx], mean, stdValue);
						std::cout << "idx\t" << idx << "\tmean\t" << mean << "\tstd\t" << stdValue << std::endl;
					}




					{
						// show distance

						std::vector<double> distances;
						for (auto& it : mergeNBNinfo.m_dis2parent) {
							if (it < 1e9) {
								distances.push_back(it);
							}
						}
						double mean(0), stdValue(0);
						calMeanAndStd(distances, mean, stdValue);
						std::cout << "distance mean\t" << mean << "\tstd\t" << stdValue << std::endl;
					}



					{

						std::vector<int> bridgeIds;
						ofec::nbn::tsp::getNarrowGap(mergeNBNinfo, bridgeIds, 0.9, 0.01);


						{
							std::ofstream out(filedir + filename + "_bridgeSolIds.txt");
							ofec::matlab::outputVector(out, bridgeIds);
							out.close();

						}
					}

				}



				//{

				//	auto solbases = nbnInfos.front().m_solbases;
				//	for (auto& it : nbnInfos.back().m_solbases) {
				//		solbases.push_back(it);
				//	}

				//	std::vector<int> filterSolId;
				//	filterUniqueSols(solbases, mergeNBNinfo.m_solbases, filterSolId, tsp_pro, rnd.get());
				//	calculateNBN(mergeNBNinfo, dataFilenames.back(), networkFilenames.back(), filename, env.get(), rnd.get());
				//}

				//std::vector<double>funnelDis;
				//std::vector<int> funnelIds;
				//calculateBelong(centerSols, mergeNBNinfo, funnelIds, funnelDis, env.get());
				//auto filedir = saveDir + "funnelInfo/";
				//std::filesystem::create_directories(filedir);

				//{
				//	ofec::ParameterVariantStream paramsStream;
				//	paramsStream << funnelIds;
				//	paramsStream << funnelDis;
				//	variants_stream::outputToFile(paramsStream, filedir + filename + "_funnelInfo.txt");
				//}


				//{
				//	std::ofstream out(filedir + filename + "funnel_network.txt");
				//	for (int idx(0); idx < funnelDis.size(); ++idx) {
				//		out << idx << "\t" << funnelIds[idx] << "\t" << funnelDis[idx] << std::endl;
				//	}
				//	out.close();
				//}

				//std::vector<std::vector<int>> belongIds(centerSols.size());
				//for (int idx(0); idx < saveDir.size(); ++idx) {
				//	if (funnelIds[idx] >= 0) {
				//		belongIds[funnelIds[idx]].push_back(idx);
				//	}
				//}

				//{
				//	{
				//		std::ofstream out(filedir + filename + "_optSolIds.txt");
				//		ofec::matlab::outputVector(out, belongIds.front());
				//		out.close();

				//	}


				//	{
				//		std::ofstream out(filedir + filename + "_furtherSolIds.txt");
				//		ofec::matlab::outputVector(out, belongIds.back());
				//		out.close();

				//	}
				//}

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


	ofec::g_working_directory = "//172.24.34.11/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "//172.24.24.151/e/DiaoYiya/code/data/ofec-data";

	//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";


	auto readdir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_nbn_rnd_lkhmore5_filterSol_final2/";
	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp2/nbn_data_eax_lkh_compare/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_subregion_two_funnel2_winServer/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_subregion_two_funnel2_linux/";
	std::filesystem::create_directories(saveDir);


	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames /*=
	{ "2202.tsp", "u574.tsp" ,"5955.tsp", "1281.tsp" ,"6702.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", }*/;

	filenames = { "5955.tsp" , "u574.tsp",   "2202.tsp", };

	//filenames = { "5955.tsp" };
	std::reverse(filenames.begin(), filenames.end());


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
			calTask(dir2, readdir2, saveDir, it);
		}
		else {
			calTask(dir1, readdir2, saveDir, it);
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