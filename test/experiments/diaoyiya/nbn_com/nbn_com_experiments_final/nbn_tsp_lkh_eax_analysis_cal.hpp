
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
#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
#include "../instance/algorithm/combination/LKH_origin/INCLUDE/LKH.h"


#include <iostream>
#include <iomanip>
#include <ctime>

#include <regex>

#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_info.h"


#include "../utility/matlab/matlab_utility.h"



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


void readTspInstance(const std::string& readDir, std::vector<std::string>& tspnames) {

	std::vector<Info> tasks;
	std::string dirpath = readDir + "/paper_com_experiment_data/totalTsp_expriment_numeric_methods/tsp_result2";


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


void ouputTask() {
	std::vector<std::string> tspnames;
	ofec::g_working_directory = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";
	readTspInstance(ofec::g_working_directory, tspnames);
	std::string saveDir = "//172.29.203.176/e/DiaoYiya/paper_com_experiment_data/experiments_added/tsp_instances/";

	{
		tspnames.resize(60);
		std::ofstream out(saveDir + "tspNames.txt");
		for (auto& it : tspnames) {
			it = it.substr(0, it.find_first_of("."));
			out << it << std::endl;
		}
		//ofec::matlab::outputVector(out, tspnames);
		out.close();
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



void calculateEAX_LKH_NBN(const std::string& readDir,
	const std::string& saveDir,
	const std::string& saveDri2,
	const std::string& tspfilename) {
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

	ofec::nbn::NBNinfo nbnInfo;
	{
		std::cout << "inputNBNinfo solutions by multithread" << std::endl;
		nbnInfo.inputNBNinfo(saveDir, tspname + "_eax_lkh");
		//nbnInfo.inputVecSolutions<int>(saveDir, tspname + "_eax_lkh", env.get());
		//nbnInfo.calculateNBNhnsw(env.get(), rnd.get());




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
			std::ofstream out2(saveDri2 + tspname + "_bestIds.txt");
			out2 << 1 << "\t" << bestIds.size() << std::endl;
			for (auto& it : bestIds) {
				out2 << it << std::endl;
			}
			out2.close();
		}


		//nbnInfo.inputTspSolutions(saveDir, tspname + "_eax_lkh", env.get());
	}


}


void calculateEAX_LKH_NBNnetwork(const std::string& readDir,
	const std::string& saveDir,
	const std::string& saveDri2,
	const std::string& tspfilename) {
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


	ofec::nbn::NBNinfo nbnInfo;
	{
		std::cout << "inputNBNinfo solutions by multithread" << std::endl;
		nbnInfo.inputNBNinfo(saveDir, tspname);
		//	nbnInfo.inputVecSolutions<int>(saveDir, tspname, env.get());
		nbnInfo.m_vFitness = UTILITY::valuesToSortedValues(nbnInfo.m_vFitness);


	}

	std::vector<std::string> filenames = { "lkh_fail", "lkh_success", "eax_fail", "eax_success" };
	std::vector<std::vector<std::vector<int>>> m_traits(4);
	std::vector<std::string> task_names = { "lkh", "eax" };
	std::vector<int> successRate;

	if (tspname == "2202") {
		std::vector<int> nearestSolIds;
		nearestSolIds.push_back(5698);

		ofec::nbn::pathToRoot(5698, nbnInfo, nearestSolIds);
		ofec::nbn::NBNinfo curInfo;
		curInfo.subset(nbnInfo, nearestSolIds);



		{
			std::ofstream out(saveDri2 + "2202" + "_eax_fail" + "_nearestPath.txt");

			outputNBNdata(out, curInfo.m_solIds, curInfo.m_belong, curInfo.m_dis2parent, curInfo.m_vFitness);;
			out.close();
		}
	}
	{
		std::vector<int> nearestSolIds;
		nearestSolIds.push_back(5698);
		std::ofstream out(saveDri2 + "2202" + "_eax_fail" + "_nearestSolIds.txt");
		ofec::matlab::outputVector(out, nearestSolIds);
		out.close();
	}



	{
		std::vector<int> nearestSolIds;
		nearestSolIds.push_back(5773);
		std::ofstream out(saveDri2 + "5955" + "_lkh_fail" + "_nearestSolIds.txt");
		ofec::matlab::outputVector(out, nearestSolIds);
		out.close();
	}



	//for (int idx(0); idx < filenames.size(); ++idx) {
	//	std::vector<int> nearestSolIds;
	//	ofec::nbn::getBestSolIds(nbnInfo, nearestSolIds);
	//	std::ofstream out(saveDri2 + tspname+ "_"+ filenames[idx] + "_nearestSolIds.txt");
	//	ofec::matlab::outputVector(out, nearestSolIds);
	//	out.close();
	//}



	/*{
		ofec::ParameterVariantStream paramsStream;
		ofec::variants_stream::inputFromFile(paramsStream, saveDir + tspname + "_trait_sucessRate.txt");
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



	}*/

	//{

	//	for (int idx(0); idx < m_traits.size(); ++idx) {
	//		std::ofstream out(saveDri2 + tspname+"_"+ filenames[idx] + "_traits.txt");
	//		auto& curTrait = m_traits[idx];
	//		std::map<int, int> mii;
	//		//	int maxIter(0);
	//		for (int itIter(0); itIter < curTrait.size(); ++itIter) {
	//			for (auto& it2 : curTrait[itIter]) {
	//				if (it2 > 0) {
	//					mii[it2] = std::max(mii[it2], itIter);
	//				}
	//			}
	//		}


	//		out << 2 << "\t" << mii.size() << std::endl;
	//		for (auto& it : mii) {
	//			out << it.first << "\t" << double(it.second) / double(curTrait.size()) << std::endl;
	//		}
	//		out.close();
	//	}




	//}

	{
		//std::vector<int> bestSolIds;
		//ofec::nbn::filterBestSolutions(bestSolIds, nbnInfo,0.99,30);
		//std::cout << "tspname\t" << tspname << "\tbest solId\t" << bestSolIds.size() << std::endl;



		//std::ofstream out(saveDri2 + tspname + "_bestSolId.txt");
		//ofec::matlab::outputVector(out, bestSolIds);
		//	out.close();
		//




	}


	{

	}




	//{

	//	int bestSolId(0);
	//	auto& fitness = nbnInfo.m_vFitness;
	//	double maxFit(0);
	//	ofec::calMax(fitness, maxFit);
	//	std::vector<int> bestIds;
	//	for (int idx(0); idx < fitness.size(); ++idx) {
	//		if (fitness[idx] == maxFit) {
	//			std::cout << idx << "\t" << maxFit << std::endl;
	//			std::cout << "maxFit\t" << maxFit << std::endl;
	//			bestIds.push_back(idx);
	//		}
	//	}
	//	auto& bestSol = nbnInfo.m_solbases[bestIds.front()];
	//	for (int idx(0); idx < 2; ++idx) {
	//		auto& curTrait = m_traits[idx];
	//		std::vector<int> lkhIds;
	//		for (auto& it : curTrait) {
	//			if (!it.empty())
	//				lkhIds.push_back(it.back());
	//		}

	//		ofec::nbn::NBNinfo curInfo;
	//		curInfo.subset(nbnInfo, lkhIds);

	//		std::cout << "update nbn NewInfo\t" << filenames[idx] << std::endl;


	//		{
	//			std::ofstream out(saveDri2 + tspname + "_" + filenames[idx] + "_lastTraitInfo.txt");

	//			outputNBNdata(out, curInfo.m_solIds, curInfo.m_belong, curInfo.m_dis2parent, curInfo.m_vFitness);;
	//			out.close();
	//		}


	//		{
	//			for (int idx(0); idx < curInfo.m_belong.size(); ++idx) {
	//				curInfo.m_dis2parent[idx] = curInfo.m_solbases[idx]->variableDistance(*bestSol, env.get());
	//				curInfo.m_belong[idx] = bestIds.front();

	//			}
	//		}
	//	}


	//	for (int idx(2); idx < 4; ++idx) {
	//		auto& curTrait = m_traits[idx];

	//		std::vector<int> eaxIds;
	//		if (!curTrait.empty()) {
	//			eaxIds = curTrait.back();
	//		}

	//		ofec::nbn::NBNinfo curInfo;
	//		curInfo.subset(nbnInfo, eaxIds);

	//		std::cout << "update nbn NewInfo\t" << filenames[idx] << std::endl;
	//		{
	//			std::ofstream out(saveDri2 + tspname + "_" + filenames[idx] + "_lastTraitInfo2.txt");

	//			//std::ofstream out(saveDri2 + tspname + "_" + filenames[idx] + "_trait_sucessRate.txt");
	//			outputNBNdata(out, curInfo.m_solIds, curInfo.m_belong, curInfo.m_dis2parent, curInfo.m_vFitness);;
	//			out.close();
	//		}
	//	}

	//}


	/*{
		int bestSolId(0);
		auto& fitness = nbnInfo.m_vFitness;
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

		auto& bestSol = nbnInfo.m_solbases[bestIds.front()];

		std::vector<std::vector<int>> nearestSolId(4);
		std::vector<double> maxDiss(4);
		std::vector<double> maxFits(4);
		for (int idx(0); idx < nearestSolId.size(); ++idx) {
			auto& curTrait = m_traits[idx];
			double maxDis = std::numeric_limits<double>::max();
			double minFit = std::numeric_limits<double>::lowest();


			nearestSolId[idx].clear();

			for (auto& it1 : curTrait) {
				for (auto& it2 : it1) {

					double curdis = bestSol->variableDistance(*nbnInfo.m_solbases[it2], env.get());
					if (curdis < maxDis) {
						maxDis = curdis;
						nearestSolId[idx].clear();
						nearestSolId[idx].push_back(it2);
					}
					else if (curdis = maxDis) {
						if (minFit < nbnInfo.m_vFitness[it2]) {
							minFit = nbnInfo.m_vFitness[it2];
							nearestSolId[idx].clear();
							nearestSolId[idx].push_back(it2);
						}
						else if (minFit == nbnInfo.m_vFitness[it2]) {
							nearestSolId[idx].push_back(it2);
						}

					}

				}
			}

			maxDiss[idx] = maxDis;
			maxFits[idx] = minFit;
		}


		for (int idx(0); idx < nearestSolId.size(); ++idx) {

			if (!nearestSolId[idx].empty()) {
				std::vector<int> path;
				ofec::nbn::pathToRoot(nearestSolId[idx].front(), nbnInfo, path);
				for (auto& it : path) {
					nearestSolId[idx].push_back(it);
				}
			}

			for (auto& it : bestIds) {
				nearestSolId[idx].push_back(it);
			}
			ofec::nbn::NBNinfo curInfo;
			curInfo.subset(nbnInfo, nearestSolId[idx]);
			std::ofstream out(saveDri2 + tspname + "_" + filenames[idx] + "_nearestBestSolInfo.txt");
			out << "update nbn NewInfo\t" << filenames[idx] << "\tMaxDis\t" << maxDiss[idx] << "\tmaxFit\t" << maxFits[idx] << std::endl;
			{


				outputNBNdata(out, curInfo.m_solIds, curInfo.m_belong, curInfo.m_dis2parent, curInfo.m_vFitness);;



			}

			out.close();
		}



	}*/










	//{
	//	std::ofstream out(saveDri2 + tspname + "_trait_sucessRate.txt");
	//	for (auto& it : task_names) {
	//		out << it << "\t";
	//	}
	//	out << std::endl;
	//	for (auto& it : successRate) {
	//		out << it << "\t";
	//	}
	//	out << std::endl;
	//	out.close();
	//}



	//{
	//	std::ofstream out(saveDri2 + tspname + "_traits.txt");
	//	for (auto& curTrait : m_traits) {
	//		std::map<int, int> mii;
	//		//	int maxIter(0);

	//		for (int itIter(0); itIter < curTrait.size(); ++itIter) {
	//			for (auto& it2 : curTrait[itIter]) {
	//				if (it2 > 0) {
	//					mii[it2] = std::max(mii[it2], itIter);
	//				}
	//			}
	//		}

	//		int maxValue(1);
	//		for (auto& it : mii) {
	//			maxValue = std::max(it.second, maxValue);
	//		}


	//		out << 2 << "\t" << mii.size() << std::endl;
	//		for (auto& it : mii) {
	//			out << it.first << "\t" << double(it.second) / double(maxValue) << std::endl;
	//		}

	//	}

	//	out.close();
	//}


	//{
	//	std::ofstream out(saveDri2 + tspname + "_nodeFit.txt");
	//	auto& fitness = nbnInfo.m_vFitness;
	//	ofec::dataNormalize(fitness);
	//	out << 1 << "\t" << fitness.size() << std::endl;
	//	for (auto& it : fitness) {
	//		out << it << std::endl;
	//	}
	//	out.close();



	//	double maxFit(0);
	//	ofec::calMax(fitness, maxFit);
	//	std::vector<int> bestIds;
	//	for (int idx(0); idx < fitness.size(); ++idx) {
	//		if (fitness[idx] == maxFit) {
	//			std::cout << idx << "\t" << maxFit << std::endl;
	//			std::cout << "maxFit\t" << maxFit << std::endl;
	//			bestIds.push_back(idx);
	//		}
	//	}
	//	{
	//		std::ofstream out2(saveDri2 + tspname + "_bestIds.txt");
	//		out2 << 1 << "\t" << bestIds.size() << std::endl;
	//		for (auto& it : bestIds) {
	//			out2 << it << std::endl;
	//		}
	//		out2.close();
	//	}
	//}





	////auto eval_fun =
	////	[](ofec::SolutionBase& sol, ofec::Environment* env) {
	////	using namespace ofec;
	////	sol.evaluate(env, false);
	////	ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
	////	sol.setFitness(pos * sol.objective(0));
	////	};

	////UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.solbases, env.get(), eval_fun);

	//std::vector<int> solIds(nbnInfo.m_belong.size());
	//for (int idx(0); idx < solIds.size(); ++idx) {
	//	solIds[idx] = idx;
	//}
	//std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	//ofec::transferNBNData(nbn_data, solIds, nbnInfo.m_belong, nbnInfo.m_dis2parent, nbnInfo.m_vFitness);

	//ofec::TreeGraphSimple nbn_graph;
	//nbn_graph.setNBNdata(nbn_data);
	////	nbn_graph.modifyBestSols(rnd.get(), env->problem(), nbnInfo.solbases);
	//nbn_graph.calNetwork(tsp_pro->numberVariables());


	//std::string savepath = saveDri2 + tspname + "_network.txt";
	//std::ofstream out(savepath);
	//nbn_graph.outputNBNnetwork(out);
	//out.close();
	//ouputNBNdata(saveDri2 + tspname + "_nbn.txt", nbn_graph.get_nbn_data());





}






void runTask_visualization() {
	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.24.34.11/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data";
	ofec::g_working_directory = "//172.24.24.151/e/DiaoYiya/code/data/ofec-data/";



	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn2/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_nbn_rnd_lkhmore5/";

	auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn_network2_1/";
	saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_nbn_rnd_lkhmore5_network_analysis/";

	std::filesystem::create_directories(saveDir);
	std::filesystem::create_directories(saveDir2);

	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";


	std::vector<std::string> filenames =
	{ "2202.tsp","1281.tsp" ,"u574.tsp" ,"6717.tsp" , "6702.tsp" , "7310.tsp", "9225.tsp", "5955.tsp", };
	//	filenames.clear();
	//	filenames = { "2202.tsp" };

	filenames.clear();

	//std::vector<Info> tasks;
	//readTspInstance(ofec::g_working_directory, filenames);
	//for (auto& it : tasks) {
	//	it.m_tspname = it.m_tspname.substr(0, it.m_tspname.find_first_of("."));
	//}



	filenames =
	{ "2202.tsp","u574.tsp" , "5955.tsp", };

	//	std::reverse(filenames.begin(), filenames.end());

	for (auto& it : filenames) {
		if (it == "u574.tsp") {
			calculateEAX_LKH_NBNnetwork(dir2, saveDir, saveDir2, it);
		}
		else {
			calculateEAX_LKH_NBNnetwork(dir1, saveDir, saveDir2, it);
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

		//runTask_visualization();

		std::cout << "runNBNvisualizationMultiTask single thread" << std::endl;
		//runNBNvisualizationMultiTask();
		runTask_visualization();
		//runTask_visualization();
//runTask();

//std::ofstream out("//172.29.203.176/e/DiaoYiya/paper_com_experiment_data/oneMaxAnalysisData/tspAlgDataAnalysis10.txt");

//std::cout << "analysis NBN tsp algorithm v1" << std::endl;
//auto readDir = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";
//analyseTSP_alg(readDir, out);
////readFiles("//172.29.203.176/e/DiaoYiya/paper_com_experiment_data/oneMax/onemaxNeighborFilterNetwork2", out);

//out.close();
	}


}