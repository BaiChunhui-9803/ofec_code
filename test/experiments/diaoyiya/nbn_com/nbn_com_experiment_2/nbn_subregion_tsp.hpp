
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

		cursol.reset(new ofec::TravellingSalesman::SolutionType(sol));
		//	cursol.reset(pro->createSolution());

	//	cursol.reset(pro->createSolution(sol));
		auto& cursolt = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*cursol);
		for (int idy(0); idy < neighbor_k; ++idy) {
			a = rnd->uniform.nextNonStd<int>(0, dim);
			b = rnd->uniform.nextNonStd<int>(0, dim);
			std::swap(cursolt.variable()[a], cursolt.variable()[b]);
		}
		cursol->evaluate(env, false);
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
	//	std::cout << "dim\t" << dim << std::endl;
	std::shared_ptr<ofec::SolutionBase> cursol;
	for (int idx(from); idx < to; ++idx) {
		cursol.reset(new ofec::TravellingSalesman::SolutionType(pro->numberObjectives(), pro->numberConstraints(), dim));
		//	cursol.reset(pro->createSolution());
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
	ofec::variants_stream::parameterStream2stringstream(paramsStream, buf);
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
	static std::string getTspSolutionFileName(const std::string& saveDir, const std::string& filename) {
		return saveDir + filename + "_solVariables.txt";
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

		std::string filepath = { NBNinfo::getTspSolutionFileName(saveDir, tspname + "_randomSamples") };
		bool filepathExists = std::filesystem::exists(filepath);
		//if (!filepathExists) {
		//	std::vector<std::shared_ptr<ofec::SolutionBase>> samples;
		//	std::cout << "generate solutions by multithread randomSample" << std::endl;
		//	sampleSolsRandomMultiThread(samples, numSamples, env.get(), rnd.get());
		//	std::vector<SolutionBase*> solbases;
		//	for (auto& it : samples) {
		//		solbases.push_back(it.get());
		//	}
		//	std::vector<int> filterSolId;
		//	NBNinfo nbnInfo;
		//	std::cout << "filterUniqueSols solutions by multithread" << std::endl;
		//	filterUniqueSols(solbases, nbnInfo.solbases, filterSolId, tsp_pro, rnd.get());
		//	UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.solbases, env.get(), eval_fun);

		//	if (nbnInfo.solbases.size() > nbnSamples)
		//		nbnInfo.solbases.resize(nbnSamples);
		//	std::cout << "calculateNBN solutions by multithread" << std::endl;
		//	calculateNBN(nbnInfo, env.get(), rnd.get());

		//	std::cout << "outputNBNinfo solutions by multithread" << std::endl;
		//	nbnInfo.outputNBNinfo(saveDir, tspname + "_randomSamples");
		//	nbnInfo.outputTspSolutions(saveDir, tspname + "_randomSamples");
		//}
	}


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

		std::vector<int> filterSolId;
		NBNinfo nbnInfo;



		filterUniqueSols(solbases, nbnInfo.solbases, filterSolId, tsp_pro, rnd.get());

		UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.solbases, env.get(), eval_fun);

		bestSol = nbnInfo.solbases.front()->getSharedPtr();
		for (auto& it : nbnInfo.solbases) {
			if (bestSol->fitness() < it->fitness()) {
				bestSol = it->getSharedPtr();
			}
		}


		std::string filepath = { saveDir + tspname + "_eaxRun_solIds.txt" };
		bool filepathExists = std::filesystem::exists(filepath);
		if (!filepathExists) {

			for (auto& it : solIds) {
				for (auto& it2 : it) {
					for (auto& it3 : it2) {
						it3 = filterSolId[it3];
					}
				}
			}

			//calculateNBN(nbnInfo, env.get(), rnd.get());
			//nbnInfo.outputNBNinfo(saveDir, tspname + "_eaxRun");
			//nbnInfo.outputTspSolutions(saveDir, tspname + "_eaxRun");
			//{

			//	ofec::ParameterVariantStream paramsStream;
			//	paramsStream << solIds.size();
			//	for (const auto& it : solIds) {
			//		paramsStream << it.size();
			//		for (const auto& it2 : it) {
			//			paramsStream << it2;
			//		}
			//	}
			//	outputToFile(paramsStream, saveDir + tspname + "_eaxRun_solIds.txt");
			//}
		}

	}


	std::vector<int> neighbors = { 100, 50,25,12,6,3 };
	auto& bestTspSol = dynamic_cast<TravellingSalesman::SolutionType&>(*bestSol);
	for (auto& neighborK : neighbors) {

		std::cout << "generate solutions by multithread sampling neighborK\t" << neighborK << std::endl;


		std::string filepath = { NBNinfo::getTspSolutionFileName(saveDir,tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK)) };
		bool filepathExists = std::filesystem::exists(filepath);

		if (!filepathExists) {
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
			UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.solbases, env.get(), eval_fun);
			calculateNBN(nbnInfo, env.get(), rnd.get());
			nbnInfo.outputNBNinfo(saveDir, tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK));
			//	nbnInfo.outputTspSolutions(saveDir, tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK));

		}
	}

}




struct Info {
	int m_suc_lkh = 0;
	int m_suc_eax = 0;
	std::string m_tspname;
	// 重载 < 运算符以提供比较方法
	bool operator<(const Info& other) const {
		// 首先比较 m_suc_lkh
		if (m_suc_lkh < other.m_suc_lkh) return true;
		if (m_suc_lkh > other.m_suc_lkh) return false;

		// m_suc_lkh 相等，则比较 m_suc_eax
		if (m_suc_eax < other.m_suc_eax) return true;
		if (m_suc_eax > other.m_suc_eax) return false;

		// m_suc_eax 也相等，则比较 m_tspname
		return m_tspname < other.m_tspname;
	}


};


void readTxtData(const std::string& filepath, std::map<std::pair<int, int>, Info>& totalInfos) {


	std::string tspname;
	int a(0), b(0);

	std::ifstream in(filepath);
	std::pair<int, int> key;
	Info curinfo;
	while (in >> tspname >> a >> b) {
		//	tspname = tspname.substr(0, tspname.find_first_of("."));
		//	std::cout << tspname << "\t" << a << "\t" << b << std::endl;
		key.first = a;
		key.second = b;
		curinfo.m_suc_lkh = a;
		curinfo.m_suc_eax = b;
		curinfo.m_tspname = tspname;
		totalInfos[key] = curinfo;
	}
	in.close();

	//int stop = -1;
}




void readTspInstance1(std::vector<std::string>& tspnames) {

	std::string dirpath = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp_expriment_numeric_methods/tsp_result2";


	std::map<std::pair<int, int>, Info> totalInfos;

	namespace fs = std::filesystem;
	// Iterate over the directory
	for (const auto& entry : fs::recursive_directory_iterator(dirpath)) {
		if (fs::is_regular_file(entry)) {
			std::cout << "File found: " << entry.path() << std::endl;

			readTxtData(entry.path().string(), totalInfos);
		}
	}

	//std::map<std::pair<int, int>, Info> totalInfos2;
	//readTxtData(dirpath + "/experiments.txt", totalInfos2);

	//for (auto& it : totalInfos2) {
	//	totalInfos[it.first] = it.second;
	//}


	//std::set<Info> sortedInfos;
	//std::vector<Info> tasks;




	//std::map<int, Info> firstInfos;
	//for (auto& it : totalInfos) {
	//	firstInfos[it.first.first] = it.second;
	//}

	//for (auto& it : firstInfos) {
	//	sortedInfos.insert(it.second);
	//	tasks.push_back(it.second);
	//	//	tspnames.push_back(it.second.m_tspname);
	//}

	//firstInfos.clear();

	//for (auto& it : totalInfos) {
	//	firstInfos[it.first.second] = it.second;
	//}

	//for (auto& it : firstInfos) {
	//	sortedInfos.insert(it.second);
	//	tasks.push_back(it.second);
	//	//	tspnames.push_back(it.second.m_tspname);
	//}


	//std::set<Info> setA;
	//std::set<Info> setB;

	//for (auto& it : totalInfos) {
	//	setA.insert(it.second);
	//}
	//for (auto& it : totalInfos2) {
	//	setB.insert(it.second);
	//}



	//std::set<Info> difference; // 用于存储差集的结果

	//// 将存在于 setA 中但不在 setB 中的所有元素插入到 difference 中
	//std::set_difference(setA.begin(), setA.end(),
	//	setB.begin(), setB.end(),
	//	std::inserter(difference, difference.begin()));


	//std::set<Info> difference2; // 用于存储差集的结果

	//// 将存在于 setA 中但不在 setB 中的所有元素插入到 difference 中
	//std::set_difference(difference.begin(), difference.end(),
	//	sortedInfos.begin(), sortedInfos.end(),
	//	std::inserter(difference2, difference2.begin()));



	//for (auto& it : difference2) {
	//	tasks.push_back(it);
	//}



	tspnames.clear();
	//std::ofstream out(dirpath + "/experiments2.txt");
	for (auto& it : totalInfos) {
		//out << it.m_tspname << "\t" << it.m_suc_lkh << "\t" << it.m_suc_eax << std::endl;
		std::cout << it.second.m_tspname << "\t" << it.second.m_suc_lkh << "\t" << it.second.m_suc_eax << std::endl;
		tspnames.push_back(it.second.m_tspname);
	}
	//out.close();




}



void readTspInstance(std::vector<std::string>& tspnames) {

	std::string dirpath = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp_expriment_numeric_methods/tsp_result2";


	std::map<std::pair<int, int>, Info> totalInfos;

	namespace fs = std::filesystem;
	// Iterate over the directory
	for (const auto& entry : fs::recursive_directory_iterator(dirpath)) {
		if (fs::is_regular_file(entry)) {
			std::cout << "File found: " << entry.path() << std::endl;

			readTxtData(entry.path().string(), totalInfos);
		}
	}

	std::map<std::pair<int, int>, Info> totalInfos2;
	readTxtData(dirpath + "/experiments.txt", totalInfos2);

	for (auto& it : totalInfos2) {
		totalInfos[it.first] = it.second;
	}


	std::set<Info> sortedInfos;
	std::vector<Info> tasks;




	std::map<int, Info> firstInfos;
	for (auto& it : totalInfos) {
		firstInfos[it.first.first] = it.second;
	}

	for (auto& it : firstInfos) {
		sortedInfos.insert(it.second);
		tasks.push_back(it.second);
		//	tspnames.push_back(it.second.m_tspname);
	}

	firstInfos.clear();

	for (auto& it : totalInfos) {
		firstInfos[it.first.second] = it.second;
	}

	for (auto& it : firstInfos) {
		sortedInfos.insert(it.second);
		tasks.push_back(it.second);
		//	tspnames.push_back(it.second.m_tspname);
	}


	std::set<Info> setA;
	std::set<Info> setB;

	for (auto& it : totalInfos) {
		setA.insert(it.second);
	}
	for (auto& it : totalInfos2) {
		setB.insert(it.second);
	}



	std::set<Info> difference; // 用于存储差集的结果

	// 将存在于 setA 中但不在 setB 中的所有元素插入到 difference 中
	std::set_difference(setA.begin(), setA.end(),
		setB.begin(), setB.end(),
		std::inserter(difference, difference.begin()));


	std::set<Info> difference2; // 用于存储差集的结果

	// 将存在于 setA 中但不在 setB 中的所有元素插入到 difference 中
	std::set_difference(difference.begin(), difference.end(),
		sortedInfos.begin(), sortedInfos.end(),
		std::inserter(difference2, difference2.begin()));



	for (auto& it : difference2) {
		tasks.push_back(it);
	}



	tspnames.clear();
	//std::ofstream out(dirpath + "/experiments2.txt");
	for (auto& it : tasks) {
		//out << it.m_tspname << "\t" << it.m_suc_lkh << "\t" << it.m_suc_eax << std::endl;
		std::cout << it.m_tspname << "\t" << it.m_suc_lkh << "\t" << it.m_suc_eax << std::endl;
		tspnames.push_back(it.m_tspname);
	}
	//out.close();




}




void runTask() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";

	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp2/nbn_data_eax_lkh_compare/";

	std::filesystem::create_directories(saveDir);


	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames /*=
	{ "2202.tsp", "u574.tsp" ,"5955.tsp", "1281.tsp" ,"6702.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", }*/;

	std::cout << "calculating nbn subresion all task ver3" << std::endl;
	readTspInstance(filenames);

	//size_t K = 7; // 假设我们想要删除前3个元素

	//if (K <= filenames.size()) {
	//	filenames.erase(filenames.begin(), std::next(filenames.begin(), K));
	//}
	// 
	// 
	filenames.resize(31);
	std::reverse(filenames.begin(), filenames.end());
	//std::reverse(filenames.begin(), filenames.end());

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
	void run(int argc, char* argv[]) {

		using namespace ofec;
		using namespace std;

		registerInstance();


		runTask();

	}


}