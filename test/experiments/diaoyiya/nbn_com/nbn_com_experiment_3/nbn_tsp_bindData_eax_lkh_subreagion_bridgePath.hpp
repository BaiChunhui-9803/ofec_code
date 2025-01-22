
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
#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"


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


void calNearestBestSolutionsThreadTask(
	int curSolId,
	int& bestSolId,
	double& bestSolDis,
	double& bestSolFit,
	const std::pair<int, int>& from_to,
	const ofec::nbn::NBNinfo& traitNBNinfo, const ofec::nbn::NBNinfo& subRegionNBNinfo, ofec::Environment* env, ofec::Random* rnd) {
	if (from_to.first >= from_to.second)return;
	auto cursol = traitNBNinfo.m_solbases[curSolId];
	bestSolId = from_to.first;
	bestSolDis = cursol->variableDistance(subRegionNBNinfo.m_solbases[from_to.first]->variableBase(), env);
	bestSolFit = subRegionNBNinfo.m_solbases[from_to.first]->fitness();
	for (int idx(from_to.first + 1); idx < from_to.second; ++idx) {
		double curDis = cursol->variableDistance(subRegionNBNinfo.m_solbases[idx]->variableBase(), env);
		if (bestSolDis > curDis) {
			bestSolDis = curDis;
			bestSolId = idx;
			bestSolFit = subRegionNBNinfo.m_solbases[idx]->fitness();
		}

		else if (bestSolDis == curDis) {
			if (bestSolFit < subRegionNBNinfo.m_solbases[idx]->fitness()) {
				bestSolDis = curDis;
				bestSolId = idx;
				bestSolFit = subRegionNBNinfo.m_solbases[idx]->fitness();
			}
			else if (bestSolFit == subRegionNBNinfo.m_solbases[idx]->fitness() && rnd->uniform.next() < 0.5) {
				bestSolDis = curDis;
				bestSolId = idx;
				bestSolFit = subRegionNBNinfo.m_solbases[idx]->fitness();
			}

			//bestSolDis = curDis;
			//bestSolId = idx;
		}
	}

}

void mergeInfo(
	int& belong,
	double& dis2best,
	double& bestFit,
	const std::vector<int>& belongs, const std::vector<double>& dis2bests, const std::vector<double>& bestFits, ofec::Random* rnd) {

	belong = belongs.front();
	dis2best = dis2bests.front();
	bestFit = bestFits.front();
	for (int idx(1); idx < belongs.size(); ++idx) {
		if (dis2best > dis2bests[idx]) {
			dis2best = dis2bests[idx];
			belong = belongs[idx];
			bestFit = bestFits[idx];
		}

		else if (dis2best == dis2bests[idx]) {
			if (bestFit < bestFits[idx]) {
				dis2best = dis2bests[idx];
				belong = belongs[idx];
				bestFit = bestFits[idx];
			}
			else if (bestFit == bestFits[idx] && rnd->uniform.next() < 0.5) {
				dis2best = dis2bests[idx];
				belong = belongs[idx];
				bestFit = bestFits[idx];
			}

			//bestSolDis = curDis;
			//bestSolId = idx;
		}
	}
}

template<typename T>
void display(const std::vector<T>& data) {
	for (auto& it : data) {
		std::cout << it << "\t";
	}
	std::cout << std::endl;
}




void getNearestBestSolutions(ofec::nbn::NBNinfo& curNBNinfo, const ofec::nbn::NBNinfo& traitNBNinfo, const ofec::nbn::NBNinfo& subRegionNBNinfo, ofec::Environment* env) {

	auto rnd = std::make_shared<ofec::Random>(0.5);


	curNBNinfo.initialize(traitNBNinfo.m_solbases.size());
	int num_task = std::thread::hardware_concurrency();
	int num_samples = subRegionNBNinfo.m_solbases.size();
	std::vector<int> belongs(num_task);
	std::vector<double> dis2bests(num_task);
	std::vector<double> bestFits(num_task);
	std::vector<std::shared_ptr<ofec::Random>> rnds(num_task);
	for (auto& it : rnds) {
		it.reset(new ofec::Random(rnd->uniform.next()));
	}

	std::vector<int> tasks;
	UTILITY::assignThreads(num_samples, num_task, tasks);
	std::pair<int, int> from_to;
	for (int curSolId(0); curSolId < traitNBNinfo.m_solbases.size(); ++curSolId) {
		std::vector<std::thread> thrds;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			thrds.push_back(std::thread(
				calNearestBestSolutionsThreadTask, curSolId, std::ref(belongs[i]), std::ref(dis2bests[i]), std::ref(bestFits[i]),
				from_to, std::cref(traitNBNinfo), std::ref(subRegionNBNinfo), env, rnds[i].get()));

		};
		for (auto& thrd : thrds)
			thrd.join();

		//std::cout << "belong" << std::endl;
		//display(belongs);
		//std::cout << "dis2bests" << std::endl;
		//display(dis2bests);
		//std::cout << "bestFits" << std::endl;
		//display(bestFits);


		mergeInfo(
			curNBNinfo.m_belong[curSolId],
			curNBNinfo.m_dis2parent[curSolId],
			curNBNinfo.m_vFitness[curSolId],
			belongs, dis2bests, bestFits, rnd.get());

		//	std::cout << "idx\t" << curSolId << "\tbelong\t" << curNBNinfo.m_belong[curSolId] << "\tm_dis2parent\t" << curNBNinfo.m_dis2parent[curSolId] << "\tm_vFitness\t" << curNBNinfo.m_vFitness[curSolId] << std::endl;



	}



}


void getBridgeData(const std::vector<int>& bestSolIds,
	const ofec::nbn::NBNinfo& subRegionNBNinfo,
	std::vector<std::vector<int>>& paths) {


	//	std::set<int> bridgeSet;


	for (auto& it : bestSolIds) {
		//bridgeSet.insert(it);
		std::vector<int> path;
		ofec::nbn::pathToRoot(it, subRegionNBNinfo, path);
		paths.push_back(path);
		//for (auto& it2 : path) {
		//	bridgeSet.insert(it2);
		//}
	}
	//
	//for (auto& it : bridgeSet) {
	//	bridgeData.push_back(it);
	//}
}


void getUniqueData(std::vector<std::vector<int>>& paths, std::vector<int>& bridgeDatas) {
	std::set<int> bridgeSet;
	for (auto& it : paths) {
		for (auto& it2 : it) {
			bridgeSet.insert(it2);
		}
	}

	for (auto& it : bridgeSet) {
		bridgeDatas.push_back(it);
	}


}

void bindDataTraitSubRegionTask(
	ofec::nbn::NBNinfo& traitNBNinfo,
	const std::string& subRegionSolDir,
	const std::string& subRegionNBNinfoDir,
	const std::string& saveDir,
	const std::string& filename,
	ofec::Environment* env
) {
	using namespace ofec;
	using namespace std;
	using namespace nbn;
	using namespace tsp;

	std::string filepath = { nbn::NBNinfo::getSolutionsFileName(subRegionSolDir,filename) };
	if (std::filesystem::exists(filepath)) {



		ofec::nbn::NBNinfo subRegionNBNinfo;
		subRegionNBNinfo.inputVecSolutions<int>(subRegionSolDir, filename, env);

		subRegionNBNinfo.inputNBNinfo(subRegionNBNinfoDir, filename);

		std::cout << "subRegionNBNinfo size" << std::endl;

		std::cout << subRegionNBNinfo.m_solbases.size() << std::endl;

		ofec::nbn::NBNinfo curNBNinfo;
		getNearestBestSolutions(curNBNinfo, traitNBNinfo, subRegionNBNinfo, env);

		std::vector<int> solIds(curNBNinfo.m_belong.size());
		for (int idx(0); idx < solIds.size(); ++idx) {
			solIds[idx] = idx;
		}
		{
			std::ofstream out(saveDir + filename + "_nearestBestSols_network.txt");
			outputNBNdata(out, solIds, curNBNinfo.m_belong, curNBNinfo.m_dis2parent, curNBNinfo.m_vFitness);
			out.close();
		}
		{
			std::ofstream out(saveDir + filename + "_nearestBestSolsIds.txt");
			ofec::matlab::outputVector(out, curNBNinfo.m_belong);
			out.close();

		}

		{
			std::vector<std::vector<int>> paths;

			getBridgeData(curNBNinfo.m_belong, subRegionNBNinfo, paths);

			NBNinfo bridgeInfo;

			{
				std::vector<int> bridgeDatas;
				getUniqueData(paths, bridgeDatas);
				bridgeInfo.subset(subRegionNBNinfo, bridgeDatas);
				std::ofstream out(saveDir + filename + "_nearestBestBridgeData.txt");
				ofec::matlab::outputVector(out, bridgeDatas);
				out.close();
			}
			{
				std::ofstream out(saveDir + filename + "_nearestBestSolPaths.txt");
				for (int idx(0); idx < paths.size(); ++idx) {
					out << idx << "\t" << curNBNinfo.m_belong[idx] << "\t" << paths[idx].size() << "\t";
					for (auto& it2 : paths[idx]) {
						out << it2 << "\t";
					}
					out << std::endl;
				}
				out.close();
			}
			{
				{
					std::ofstream out(saveDir + filename + "_bridgeSols_network.txt");
					outputNBNdata(out, solIds, bridgeInfo.m_belong, bridgeInfo.m_dis2parent, bridgeInfo.m_vFitness);
					out.close();
				}
				{
					std::ofstream out(saveDir + filename + "_bridgeSolsIds.txt");
					ofec::matlab::outputVector(out, bridgeInfo.m_belong);
					out.close();

				}
			}

		}
	}
}


void bindDataTraitSubRegion(
	const std::string& readDir,
	const std::string& traitDir,
	const std::string& subRegionSolDir,
	const std::string& subRegionNBNinfoDir,
	const std::string& saveDir,
	const std::string& tspfilename
) {
	using namespace ofec;
	using namespace std;
	using namespace nbn;
	using namespace tsp;

	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));

	std::cout << "calculating filename\t" << tspname << std::endl;

	std::shared_ptr<Environment> env;
	genenrateTSPenv(readDir, tspfilename, env);


	std::string filepath;
	filepath = NBNinfo::getSolutionsFileName(traitDir, tspname + "_optSols");
	if (!std::filesystem::exists(filepath)) {
		return;
	}
	ofec::nbn::NBNinfo traitNBNinfo;
	traitNBNinfo.inputVecSolutions<int>(traitDir, tspname + "_optSols", env.get());
	std::cout << "traitNBNinfo size" << std::endl;
	std::cout << traitNBNinfo.m_solbases.size() << std::endl;
	std::string filename;
	filename = tspname + "_randomSamples";

	bindDataTraitSubRegionTask(
		traitNBNinfo, subRegionSolDir, subRegionNBNinfoDir, saveDir, filename, env.get());

	std::vector<int> neighbors = { 100, 50,25,12,6,3 };
	//std::reverse(neighbors.begin(), neighbors.end());

	for (auto& neighborK : neighbors) {

		std::cout << "get solutions by multithread sampling neighborK\t" << neighborK << std::endl;
		filename = tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK);
		bindDataTraitSubRegionTask(
			traitNBNinfo, subRegionSolDir, subRegionNBNinfoDir, saveDir, filename, env.get());



		//std::string filepath = { NBNinfo::getSolutionsFileName(subRegionDir,filename) };
		//bool filepathExists = std::filesystem::exists(filepath);
		//if (filepathExists) {
		//	ofec::nbn::NBNinfo nbnInfo;
		//	nbnInfo.inputVecSolutions<int>(subRegionDir, tspname + "_randomSamples" + "_neighborK_" + std::to_string(neighborK), env.get());
		//}
	}


}








void runTask() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
	ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";




	const std::string traitDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_nbn_hnswEqualBetter_analysis/";
	const std::string& subRegionSolDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread/";
	const std::string& subRegionNbnInfoDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread_equalBetter/";

	const std::string& saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread_bindData/";
	std::filesystem::create_directories(saveDir);


	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames /*=
	{ "2202.tsp", "u574.tsp" ,"5955.tsp", "1281.tsp" ,"6702.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", }*/;

	filenames = { "u574.tsp",  "5955.tsp" ,   "2202.tsp", };


	//std::reverse(filenames.begin(), filenames.end());


	std::cout << "calculating nbn subresion all task nbn_hnswEqualBetter_bindData" << std::endl;
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
			bindDataTraitSubRegion(dir2, traitDir, subRegionSolDir, subRegionNbnInfoDir, saveDir, it);
		}
		else {
			bindDataTraitSubRegion(dir1, traitDir, subRegionSolDir, subRegionNbnInfoDir, saveDir, it);
			//	calTask(dir1, saveDir, saveDir2, it);
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