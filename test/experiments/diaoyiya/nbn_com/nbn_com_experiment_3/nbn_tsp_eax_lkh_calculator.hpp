
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


struct LKH_RunInfo {
	std::vector<unsigned> m_seeds;
	std::vector<std::vector<std::vector<int>>> m_sols;
	std::vector<double> m_costs;
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
		std::vector<double> costs;
		lkh_alg.run(readDir, tspfilename, totalSols);

		//bool correctRecord = true;


		//for (auto& it : curInfo.m_costs) {
		//	if (lkh_alg.BestCost < it) {
		//		correctRecord = false;
		//	}
		//}

		//if (!correctRecord) {
		//	std::cout << "error at" << tspfilename << std::endl;
		//}

		for (auto& it : totalSols) {
			for (auto& it2 : it) {
				--it2;
			}
		}
	}




}


void runLKHmultiThread(
	const std::string& readDir,
	const std::string& tspfilename,
	LKH_RunInfo& totalInfo) {
	using namespace ofec;
	int numRuns(30);
	//numRuns = 1;
	totalInfo.resetSize(numRuns);

	unsigned initSeed = 1;
	unsigned SeedStep = 1e4;

	ofec::Random rnd(0.5);
	for (auto& it : totalInfo.m_seeds) {
		it = initSeed;
		initSeed += SeedStep;
		it = rnd.uniform.nextNonStd<unsigned>(1, std::numeric_limits<unsigned>::max());
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


struct TaskInfo {
	std::ofstream m_out;
	std::set<int> m_lkh_suc;
	std::set<int> m_eax_suc;

};

void runEAX_LKH(const std::string& readDir,
	const std::string& saveDir,
	const std::string& tspfilename, TaskInfo& globalInfo) {
	using namespace ofec;
	using namespace std;
	using namespace chrono;
	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));
	std::cout << "calculating filename\t" << tspname << std::endl;
	std::shared_ptr<Environment> env;
	genenrateTSPenv(readDir, tspfilename, env);

	std::string proname = "TSP";
	auto pro = env->problem();
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
	auto rnd = std::make_shared<ofec::Random>(0.5);




	LKH_RunInfo lkhInfo;
	runLKHmultiThread(readDir, tspfilename, lkhInfo);
	std::vector<std::vector<std::vector<std::shared_ptr<ofec::SolutionBase>>>> sols;
	runEax(readDir, tspfilename, sols);

	std::vector<std::vector<int>> lkh_ids;
	std::vector<std::vector<std::vector<int>>> eax_ids;

	std::vector<SolutionBase*> totalSols;
	//	int solId(0);
	eax_ids.resize(sols.size());
	for (int idx(0); idx < eax_ids.size(); ++idx) {
		eax_ids[idx].resize(sols[idx].size());
		for (int idy(0); idy < eax_ids[idx].size(); ++idy) {
			eax_ids[idx][idy].resize(sols[idx][idy].size());
			for (int idz(0); idz < eax_ids[idx][idy].size(); ++idz) {
				eax_ids[idx][idy][idz] = totalSols.size();
				totalSols.push_back(sols[idx][idy][idz].get());

			}
		}
	}

	std::vector<std::shared_ptr<SolutionBase>> lkh_sols;
	auto& lkh_solx = lkhInfo.m_sols;
	lkh_ids.resize(lkh_solx.size());
	for (int idx(0); idx < lkh_solx.size(); ++idx) {
		lkh_ids[idx].resize(lkh_solx[idx].size());
		for (int idy(0); idy < lkh_solx[idx].size(); ++idy) {
			std::shared_ptr<SolutionBase> cursol;
			cursol.reset(pro->createSolution());
			auto& tspSol = dynamic_cast<TravellingSalesman::SolutionType&>(*cursol);
			tspSol.variable().vect() = lkh_solx[idx][idy];
			lkh_sols.push_back(cursol);
			lkh_ids[idx][idy] = totalSols.size();
			totalSols.push_back(cursol.get());

		}
	}

	std::vector<SolutionBase*> filterSols;
	std::vector<int> solId2FilterId;

	nbn::tsp::filterUniqueSols(totalSols, filterSols, solId2FilterId, tsp_pro, rnd.get());

	//{
	//	// for test 
	//	std::vector<int> dupulateSol = { 0,201,78,351,237,328,382,47,151,164,73,236,277,19,463,388,2,392,128,25,184,39,360,314,207,104,175,224,409,155,343,247,13,444,61,455,119,192,103,157,244,345,329,394,74,441,332,117,427,486,307,66,208,251,308,121,152,453,58,352,331,403,418,261,83,142,231,59,233,53,304,448,173,446,35,381,225,340,216,305,395,130,69,295,145,255,219,14,359,312,46,22,206,125,268,315,43,349,327,379,457,220,398,434,187,378,302,232,355,143,293,259,383,411,167,163,71,44,217,40,487,266,112,479,484,262,489,422,459,412,414,467,186,290,384,6,488,171,67,33,483,485,102,24,211,20,118,320,473,183,153,499,292,79,282,319,363,361,87,406,393,276,325,139,471,426,288,110,26,405,154,49,273,63,174,404,182,51,300,138,64,80,86,420,365,158,107,30,323,374,15,429,126,194,94,136,373,478,474,48,75,204,460,270,433,299,494,430,29,54,492,461,245,431,36,447,341,468,5,168,432,146,188,357,12,271,185,337,339,191,98,137,390,32,166,371,482,367,177,196,131,249,195,135,132,99,200,144,239,165,34,490,283,324,149,316,93,101,439,193,18,178,72,45,23,402,368,333,55,475,37,294,437,296,134,370,115,238,197,114,230,498,120,452,242,442,129,97,369,346,289,301,122,416,425,209,272,281,326,42,159,347,56,91,417,68,123,397,21,140,419,318,52,124,190,497,298,95,172,458,28,330,274,385,203,280,250,303,170,348,234,223,213,212,7,228,92,313,464,285,362,88,141,96,462,16,321,443,336,396,106,334,309,116,465,10,450,291,60,297,156,317,57,410,375,160,258,229,9,372,189,401,108,243,286,150,495,377,380,226,77,335,267,496,284,338,493,399,148,50,11,235,180,214,263,366,344,240,113,451,353,3,445,421,469,41,89,306,350,287,82,269,81,256,248,222,408,310,481,407,221,436,424,176,449,27,376,356,17,400,127,440,31,90,456,428,472,179,169,438,342,253,480,470,435,8,415,1,241,70,205,227,466,476,210,62,264,311,100,358,279,198,76,260,85,147,322,161,477,215,389,413,202,423,111,387,278,162,65,4,199,257,109,181,38,105,391,364,246,454,84,265,133,252,386,254,491,275,218,354 };

	//	std::shared_ptr<ofec::SolutionBase> cursol(tsp_pro->createSolution());
	//	auto& tspSol = dynamic_cast<TravellingSalesman::SolutionType&>(*cursol);
	//	tspSol.variable().vect() = dupulateSol;
	//	
	//	int numdupulateSol = 0;
	//	for (auto& it : filterSols) {
	//		if (it->variableDistance(tspSol.variableBase(), env.get())==0) {
	//			++numdupulateSol;
	//		}
	//	}

	//	std::cout << "numdupulateSol\t" << numdupulateSol << std::endl;
	//}



	auto eval_fun =
		[](ofec::SolutionBase& sol, ofec::Environment* env) {
		using namespace ofec;
		sol.evaluate(env, false);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};

	UTILITY::evaluateRandomSolutionsMultiThreads(filterSols, env.get(), eval_fun);


	int solSize = filterSols.size();
	int sqrtSize = sqrt(solSize);
	solSize = sqrtSize * sqrtSize;
	int removeSize = filterSols.size() - solSize;
	removeSize = 0;
	std::vector<int> sortIds(filterSols.size());
	for (int idx(0); idx < sortIds.size(); ++idx) {
		sortIds[idx] = idx;
	}

	std::sort(sortIds.begin(), sortIds.end(), [&](int a, int b) {
		return filterSols[a]->fitness() < filterSols[b]->fitness();
		});

	sortIds.resize(removeSize);
	std::vector<bool> leftSolFlag(filterSols.size(), true);
	for (auto& it : sortIds) leftSolFlag[it] = false;
	std::vector<int> filterIdToRemoveId(filterSols.size(), -1);

	std::vector<SolutionBase*> leftSols;

	for (int idx(0); idx < filterSols.size(); ++idx) {
		if (leftSolFlag[idx]) {
			filterIdToRemoveId[idx] = leftSols.size();
			leftSols.push_back(filterSols[idx]);

		}
	}


	{
		std::vector<std::string> filenames = { "lkh_fail", "lkh_success", "eax_fail", "eax_success" };
		std::vector<std::vector<std::vector<int>>> m_traits(4);
		std::vector<std::string> task_names = { "lkh", "eax" };
		std::vector<int> successRate(2, 0);
		std::vector<bool> traitFlag(4, false);
		double bestFit(std::numeric_limits<double>::lowest());
		for (auto& it : leftSols) {
			bestFit = std::max(it->fitness(), bestFit);
		}
		for (auto& it : lkh_ids) {
			for (auto& it2 : it) {
				it2 = solId2FilterId[it2];
				it2 = filterIdToRemoveId[it2];
			}
		}

		for (auto& it : eax_ids) {
			for (auto& it2 : it) {
				for (auto& it3 : it2) {
					it3 = solId2FilterId[it3];
					it3 = filterIdToRemoveId[it3];
				}
			}
		}

		int traitId(0);
		for (auto& it : lkh_ids) {
			double bestFitRun(std::numeric_limits<double>::lowest());
			for (auto& it2 : it) {
				if (it2 >= 0) {
					bestFitRun = std::max(leftSols[it2]->fitness(), bestFitRun);
				}
			}

			traitId = -1;
			if (bestFitRun == bestFit) {
				++successRate.front();
				traitId = 1;
				if (!traitFlag[traitId]) {
					traitFlag[traitId] = true;
				}
				else traitId = -1;
			}
			else {
				traitId = 0;
				if (!traitFlag[traitId]) {
					traitFlag[traitId] = true;
				}
				else traitId = -1;
			}

			if (traitId >= 0) {
				auto& curTrait = m_traits[traitId];
				curTrait.resize(it.size());
				for (int idx(0); idx < curTrait.size(); ++idx) {
					curTrait[idx].push_back(it[idx]);
				}
				//for (auto& it2 : it) {
				//	if (it2 >= 0) curTrait.push_back(it2);
				//}
			}
		}

		for (auto& it : eax_ids) {

			double bestFitRun(std::numeric_limits<double>::lowest());
			for (auto& it2 : it) {
				for (auto& it3 : it2) {
					if (it3 >= 0) {
						bestFitRun = std::max(leftSols[it3]->fitness(), bestFitRun);
					}
				}
			}


			traitId = -1;
			if (bestFitRun == bestFit) {
				++successRate.back();
				traitId = 3;
				if (!traitFlag[traitId]) {
					traitFlag[traitId] = true;
				}
				else traitId = -1;
			}
			else {
				traitId = 2;
				if (!traitFlag[traitId]) {
					traitFlag[traitId] = true;
				}
				else traitId = -1;
			}

			if (traitId >= 0) {
				auto& curTrait = m_traits[traitId];
				curTrait = it;
				/*			std::set<int> solIdSet;
							for (auto& it2 : it) {
								for (auto& it3 : it2) {
									if (it3 >= 0)solIdSet.insert(it3);
								}
							}

							for (auto& it : solIdSet) {
								curTrait.push_back(it);
							}*/
			}

		}


		bool isOutput(false);

		{
			if (globalInfo.m_lkh_suc.find(successRate.front()) == globalInfo.m_lkh_suc.end()) {
				globalInfo.m_lkh_suc.insert(successRate.front());
				isOutput = true;
			}
			if (globalInfo.m_eax_suc.find(successRate.back()) == globalInfo.m_eax_suc.end()) {
				globalInfo.m_eax_suc.insert(successRate.back());
				isOutput = true;
			}

		}

		/*if (isOutput) */ {


			{
				globalInfo.m_out << tspfilename << "\t" << successRate.front() << "\t" << successRate.back() << std::endl;
				std::cout << tspfilename << "\t" << successRate.front() << "\t" << successRate.back() << std::endl;
			}

			{
				ofec::ParameterVariantStream paramsStream;
				paramsStream << m_traits.size();
				for (const auto& it : m_traits) {
					paramsStream << it.size();
					for (const auto& it2 : it) {
						paramsStream << it2;
					}
				}
				paramsStream << successRate;

				variants_stream::outputToFile(paramsStream, saveDir + tspname + "_trait_sucessRate.txt");
			}


			{

				nbn::NBNinfo nbnInfo;
				nbnInfo.m_solbases = leftSols;


				calculateNBN(nbnInfo, saveDir, saveDir, tspname, env.get(), rnd.get());

				//std::cout << "calculateNBN solutions by multithread" << std::endl;
				//nbnInfo.calculateNBNhnswEqualBetter(env.get(), rnd.get());

				//std::cout << "outputNBNinfo solutions by multithread" << std::endl;
				//nbnInfo.outputNBNinfo(saveDir, tspname + "_eax_lkh");
				//nbnInfo.outputVecSolutions<int>(saveDir, tspname + "_eax_lkh", env.get());
				{
					std::ofstream out(saveDir + tspname + "_traits.txt");
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
					std::ofstream out(saveDir + tspname + "_nodeFit.txt");
					auto& fitness = nbnInfo.m_vFitness;
					ofec::dataNormalize(fitness);
					out << 1 << "\t" << fitness.size() << std::endl;
					for (auto& it : fitness) {
						out << it << std::endl;
					}
					out.close();



					double maxFit(0);
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
						std::ofstream out2(saveDir + tspname + "_bestIds.txt");
						out2 << 1 << "\t" << bestIds.size() << std::endl;
						for (auto& it : bestIds) {
							out2 << it << std::endl;
						}
						out2.close();
					}
				}
			}




		}

	}









	//{
	//	ofec::EAX_TSP_trait trait;
	//	std::vector<int> cursolx;
	//	for (int idx(0); idx < leftSols.size(); ++idx) {
	//		ofec::TravellingSalesman::getSolutionX(*leftSols[idx], cursolx);
	//		//auto& tspsol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*leftSols[idx]);
	//		trait.m_sols.push_back(cursolx);
	//	}

	//	trait.outputIDEESols(saveDir + tspname + "_IDEE_solutions.txt");
	//}


}


void runTask() {
	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "//172.29.41.69/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.24.24.151/e/DiaoYiya/code/data/ofec-data";
	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_nbn_rnd_lkhmore/";

	std::filesystem::create_directories(saveDir);


	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";
	std::string dir3 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin2/";
	std::vector<std::string> filenames;
	//{ "2202.tsp","1281.tsp" ,"u574.tsp" ,"6717.tsp" , "6702.tsp" , "7310.tsp", "9225.tsp", "5955.tsp", };


	filenames = { "2202.tsp","u574.tsp" , "5955.tsp", };


	//readTspInstance2(ofec::g_working_directory, filenames);


	//readDirFun(dir1, filenames);

	//std::reverse(filenames.begin(), filenames.end());


	std::cout << "run EAX vs LKH run compareResult hnsw equal better" << std::endl;

	//filenames.resize(60);
	//filenames.erase(filenames.begin());
	//	filenames.clear();
	//	filenames = { "2202.tsp" };

		//	std::reverse(filenames.begin(), filenames.end());
	//size_t K = 11;

	//if (K <= filenames.size()) {
	//	filenames.erase(filenames.begin(), std::next(filenames.begin(), K));
	//}

	//std::reverse(filenames.begin(), filenames.end());


	TaskInfo globalInfo;
	globalInfo.m_out.open(saveDir + "tsp_instance2.txt");

	for (auto& it : filenames) {
		if (it == "u574.tsp") {
			runEAX_LKH(dir2, saveDir, it, globalInfo);
		}
		else {
			runEAX_LKH(dir1, saveDir, it, globalInfo);
		}
	}


	globalInfo.m_out.close();



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