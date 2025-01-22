
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



void analyseEAX_LKH_NBN(const std::string& readDir,
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
		nbnInfo.inputVecSolutions<int>(saveDir, tspname, env.get());
		//	nbnInfo.calculateNBNhnsw(env.get(), rnd.get());

		nbnInfo.m_vFitness = UTILITY::valuesToSortedValues(nbnInfo.m_vFitness);

		//nbnInfo.inputTspSolutions(saveDir, tspname + "_eax_lkh", env.get());
	}

	std::vector<int> bestIds;
	{
		int bestSolId(0);
		auto& fitness = nbnInfo.m_vFitness;
		double maxFit(0);
		ofec::calMax(fitness, maxFit);

		for (int idx(0); idx < fitness.size(); ++idx) {
			if (fitness[idx] == maxFit) {
				std::cout << idx << "\t" << maxFit << std::endl;
				std::cout << "maxFit\t" << maxFit << std::endl;
				bestIds.push_back(idx);
			}
		}
	}




	{

		ofec::dataNormalize(nbnInfo.m_vFitness);

		double fitnessThread = 0.99;
		double disThread = 10;


		std::vector<int> optIds;
		//std::vector<int> solId;
		//std::vector<int> belong;
		//std::vector<double> dis2parent;
		//std::vector<double> fitness;


		for (int idx(0); idx < nbnInfo.m_vFitness.size(); ++idx) {
			if (nbnInfo.m_vFitness[idx] >= fitnessThread && nbnInfo.m_dis2parent[idx] >= disThread) {
				optIds.push_back(idx);
				//solId.push_back(idx);
				//belong.push_back(nbnInfo.m_belong[idx]);
				//dis2parent.push_back(nbnInfo.m_dis2parent[idx]);
				//fitness.push_back(nbnInfo.m_vFitness[idx]);
			}
		}

		ofec::nbn::NBNinfo optNBN;
		optNBN.subset(nbnInfo, optIds);
		std::ofstream out(saveDri2 + tspname + "_optSols_network.txt");
		outputNBNdata(out, optNBN.m_solIds, optNBN.m_belong, optNBN.m_dis2parent, optNBN.m_vFitness);
		out.close();
		{
			std::ofstream out(saveDri2 + tspname + "_optSolIds.txt");
			ofec::matlab::outputVector(out, optIds);
			out.close();

		}

		//int bestSolId(0);
		//{
		//	double maxFit = std::numeric_limits<double>::lowest();
		//	for (int idx(0); idx < optNBN.m_vFitness.size(); ++idx) {
		//		if (maxFit < optNBN.m_vFitness[idx]) {
		//			maxFit = optNBN.m_vFitness[idx];
		//			bestSolId = idx;
		//		}
		//		
		//	}
		//}

		{

			auto bestSol = nbnInfo.m_solbases[bestIds.front()];
			std::vector<double> dis2opt(optNBN.m_dis2parent.size());
			std::vector<int> filterSolIds;

			for (int idx(0); idx < dis2opt.size(); ++idx) {
				dis2opt[idx] = optNBN.m_solbases[idx]->variableDistance(*bestSol, env.get());

			}

			std::vector<int> filterIds;
			for (int idx(0); idx < dis2opt.size(); ++idx) {
				if (dis2opt[idx] >= 10 && dis2opt[idx] <= 17) {
					filterIds.push_back(idx);
				}
	/*			else if (dis2opt[idx] == 0) {
					filterIds.push_back(idx);
				}*/
			}

			{
				std::vector<int> solIds;
				ofec::nbn::getBestSolIds(optNBN, solIds);
				auto bestSolId = solIds.front();
				filterIds.push_back(bestSolId);
				//centerSols.push_back(optNBNinfo.m_solbases[solIds.front()]);
			}


			ofec::nbn::NBNinfo filterNBN;
			filterNBN.subset(optNBN, filterIds);


			filterNBN.outputNBNinfo(saveDri2, tspname + "_filterOptSols");
			filterNBN.outputVecSolutions<int>(saveDri2, tspname + "_filterOptSols", env.get());
			std::ofstream out(saveDri2 + tspname + "_filterOptSols_network.txt");
			outputNBNdata(out, filterNBN.m_solIds, filterNBN.m_belong, filterNBN.m_dis2parent, filterNBN.m_vFitness);
			out.close();
			{
				std::ofstream out(saveDri2 + tspname + "_filterOptSols_Dis.txt");
				out << "CurId\tOrigin\tBelong\tDis2parent\tFitness\tDis2Opt" << std::endl;
				for (int idx(0); idx < optNBN.m_vFitness.size(); ++idx) {
					out << idx << "\t" << optNBN.m_solIds[idx] << "\t" << optNBN.m_belong[idx] << "\t" << optNBN.m_dis2parent[idx] << "\t" << optNBN.m_vFitness[idx] << "\t" << dis2opt[idx] << std::endl;
				}
				out.close();
			}
		}


		{

			optNBN.outputNBNinfo(saveDri2, tspname + "_optSols");
			optNBN.outputVecSolutions<int>(saveDri2, tspname + "_optSols", env.get());

		}
	}

}





void runTask_visualization() {
	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.24.34.11/share/2018/diaoyiya/ofec-data";
	//	ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data";
	ofec::g_working_directory = "//172.24.24.151/e/DiaoYiya/code/data/ofec-data";


	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn2/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_nbn_rnd_lkhmore5/";

	auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn_network2_1/";
	saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_nbn_rnd_lkhmore5_filterSol_final/";

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
			analyseEAX_LKH_NBN(dir2, saveDir, saveDir2, it);
		}
		else {
			analyseEAX_LKH_NBN(dir1, saveDir, saveDir2, it);
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