
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

#include <regex>

#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_info.h"


#include "../utility/matlab/matlab_utility.h"

// oneMaxInfo

struct NBNanalyseInfo {
	int EpistasisBlockSize = 0;
	int NeutralityMu = 0;
	int RuggednessGamma = 0;
	int neighborK = 0;


	ofec::nbn::NBNinfo nbnInfo;

	std::vector<double> values;




};


void analyseNameOneMax(std::vector<std::string>& values) {

	values.push_back("NBDstd");
	values.push_back("MaxNeutralBasins");

	values.push_back("numberOpts");

	values.push_back("MaxNBDofOptimum");
	values.push_back("basinSizesSTD");

	values.push_back("maxNeutralitySons");

	values.push_back("bestNBDstd");
}


void analyseOneMax(NBNanalyseInfo& info) {


	info.nbnInfo.modify();

	using namespace ofec::nbn;
	auto& values = info.values;
	auto& nbnInfo = info.nbnInfo;

	std::vector<int> neurtual_size;

	values.push_back(calculateNBDstd(nbnInfo));
	values.push_back(calculateMaxNeutralBasins(neurtual_size, nbnInfo));

	std::vector<int> optIds;
	getOptimumIds(optIds, nbnInfo, 0.02, 6);

	values.push_back(optIds.size());

	values.push_back(calculateMaxNBDofOptimum(nbnInfo, optIds));

	std::vector<std::vector<int>> basins;
	calculateBasins(basins, nbnInfo, optIds);

	values.push_back(basinSizesSTD(basins));

	values.push_back(maxNeutralitySons(nbnInfo));
	values.push_back(getNBDvalue(nbnInfo));
}




void readFiles(const std::string& directoryPath, std::ofstream& out) {

	namespace fs = std::filesystem;
	const std::string extension = "_nbn.txt";

	// 正则表达式匹配文件名后缀
	std::regex fileNameRegex(".*" + extension);
	// 用于存储找到的文件路径
	std::vector<fs::path> matchedFiles;

	// 遍历目录
	for (const auto& entry : fs::recursive_directory_iterator(directoryPath)) {
		if (fs::is_regular_file(entry)) {
			// 检查文件名是否以特定后缀结束
			if (std::regex_match(entry.path().filename().string(), fileNameRegex)) {
				matchedFiles.push_back(entry.path());
				//	std::cout << entry.path().filename()<< std::endl;
			}
		}
	}


	std::vector<std::string> featureNames = { "EpistasisBlockSize",
	"NeutralityMu","RuggednessGamma" };

	std::map<std::string, int> featureNameIds;
	for (int idx(0); idx < featureNames.size(); ++idx) {
		featureNameIds[featureNames[idx]] = idx;
	}
	std::string neighborKname = "neighborK";
	std::map<std::string, int> tags;
	std::map< std::pair<std::string, std::pair<int, int>>, NBNanalyseInfo> totalInfos;
	std::pair<std::string, std::pair<int, int>> infoKey;


	for (auto& it : matchedFiles) {
		auto name = it.filename().string();
		auto keys = UTILITY::split(name, "_");
		keys.pop_back();
		tags.clear();
		for (int idx(1); idx < keys.size(); idx += 2) {
			tags[keys[idx]] = std::stoi(keys[idx + 1]);
		}
		int keyFeature = -1;
		for (int idx(0); idx < featureNames.size(); ++idx) {
			if (tags[featureNames[idx]] != 0) {
				keyFeature = idx;
				break;
			}
		}


		std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
		ofec::inputNBNdata(it.string(), nbn_data);
		NBNanalyseInfo oneMaxInfo;
		oneMaxInfo.nbnInfo.fromNBNdata(nbn_data);
		analyseOneMax(oneMaxInfo);

		if (keyFeature == -1) {


			for (keyFeature = 0; keyFeature < featureNames.size(); ++keyFeature) {
				infoKey.first = featureNames[keyFeature];
				infoKey.second.first = tags[infoKey.first];
				infoKey.second.second = tags[neighborKname];
				totalInfos[infoKey] = oneMaxInfo;
			}
		}
		else {
			infoKey.first = featureNames[keyFeature];
			infoKey.second.first = tags[infoKey.first];
			infoKey.second.second = tags[neighborKname];
			totalInfos[infoKey] = oneMaxInfo;
		}

		//int stop = -1;

	}


	std::vector<std::set<int>> featureRanges(featureNames.size());

	std::set<int> neighborRanges;

	for (auto& it : totalInfos) {
		featureRanges[featureNameIds[it.first.first]].insert(it.first.second.first);
		neighborRanges.insert(it.first.second.second);


	}



	std::vector<std::string> names;
	analyseNameOneMax(names);
	for (int idx(0); idx < featureNames.size(); ++idx) {
		out << featureNames[idx] << "\t";
		for (auto& neighborK : neighborRanges) {
			out << neighborKname << "\t";
			for (auto& it : names) {
				out << it << "\t";
			}
		}

		out << std::endl;


		infoKey.first = featureNames[idx];
		for (auto& featureValue : featureRanges[idx]) {
			out << featureValue << "\t";

			infoKey.second.first = featureValue;
			for (auto& neighborK : neighborRanges) {
				out << neighborK << "\t";
				infoKey.second.second = neighborK;

				auto& oneMaxInfo = totalInfos[infoKey];
				for (auto& it : oneMaxInfo.values) {
					out << it << "\t";
				}
			}
			out << std::endl;

		}


	}

	//int stop = -1;



}


struct Info {

	std::string m_dir;


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


void readTxtData(const std::string& filepath, std::vector<Info>& totalInfos) {


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
		totalInfos.push_back(curinfo);
		//totalInfos[key] = curinfo;
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


//void analyseTSP_alg(const std::string& readDir, std::ofstream& out) {
//	std::vector<Info> tasks;
//	readTspInstance(readDir, tasks);
//	for (auto& it : tasks) {
//		it.m_tspname = it.m_tspname.substr(0, it.m_tspname.find_first_of("."));
//	}
//
//
//	std::vector<NBNanalyseInfo> nbnAnInfo(tasks.size());
//
//	std::string nbnReadFile = readDir + +"/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn/";
//
//	for (int idx(0); idx < nbnAnInfo.size(); ++idx) {
//
//		nbnAnInfo[idx].nbnInfo.inputNBNinfo(nbnReadFile, tasks[idx].m_tspname + "_eax_lkh");
//		//std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
//		//ofec::inputNBNdata(tasks[idx].m_tspname+"", nbn_data);
//		//NBNanalyseInfo oneMaxInfo;
//		//oneMaxInfo.nbnInfo.fromNBNdata(nbn_data);
//		analyseOneMax(nbnAnInfo[idx]);
//
//	}
//
//
//	std::vector<std::string> headNames = { "tspname" ,"suc_lkh" ,"suc_eax" };
//	std::vector<std::string> headNames2;
//
//	analyseNameOneMax(headNames2);
//	for (auto& it : headNames) {
//		out << it << "\t";
//	}
//	for (auto& it : headNames2) {
//		out << it << "\t";
//	}
//	out << std::endl;
//
//	for (int idx(0); idx < nbnAnInfo.size(); ++idx) {
//		auto& it = tasks[idx];
//		out << it.m_tspname << "\t" << it.m_suc_lkh << "\t" << it.m_suc_eax << "\t";
//		auto& values = nbnAnInfo[idx].values;
//		for (auto& it2 : values) {
//			out << it2 << "\t";
//		}
//		out << std::endl;
//	}
//}


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
		nbnInfo.inputNBNinfo(saveDir, tspname + "_eax_lkh");
		//nbnInfo.inputTspSolutions(saveDir, tspname + "_eax_lkh", env.get());
	}

	std::vector<std::string> filenames = { "lkh_fail", "lkh_success", "eax_fail", "eax_success" };
	std::vector<std::vector<std::vector<int>>> m_traits(4);
	std::vector<std::string> task_names = { "lkh", "eax" };
	std::vector<int> successRate;

	{
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


	}



	{
		std::ofstream out(saveDri2 + tspname + "_trait_sucessRate.txt");
		for (auto& it : task_names) {
			out << it << "\t";
		}
		out << std::endl;
		for (auto& it : successRate) {
			out << it << "\t";
		}
		out << std::endl;
		out.close();
	}



	{
		std::ofstream out(saveDri2 + tspname + "_traits.txt");
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
		std::ofstream out(saveDri2 + tspname + "_nodeFit.txt");
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
			std::ofstream out2(saveDri2 + tspname + "_bestIds.txt");
			out2 << 1 << "\t" << bestIds.size() << std::endl;
			for (auto& it : bestIds) {
				out2 << it << std::endl;
			}
			out2.close();
		}
	}





	//auto eval_fun =
	//	[](ofec::SolutionBase& sol, ofec::Environment* env) {
	//	using namespace ofec;
	//	sol.evaluate(env, false);
	//	ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
	//	sol.setFitness(pos * sol.objective(0));
	//	};

	//UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.solbases, env.get(), eval_fun);

	std::vector<int> solIds(nbnInfo.m_belong.size());
	for (int idx(0); idx < solIds.size(); ++idx) {
		solIds[idx] = idx;
	}
	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	ofec::transferNBNData(nbn_data, solIds, nbnInfo.m_belong, nbnInfo.m_dis2parent, nbnInfo.m_vFitness);

	ofec::TreeGraphSimple nbn_graph;
	nbn_graph.setNBNdata(nbn_data);
	//	nbn_graph.modifyBestSols(rnd.get(), env->problem(), nbnInfo.solbases);
	nbn_graph.calNetwork(tsp_pro->numberVariables());


	std::string savepath = saveDri2 + tspname + "_network.txt";
	std::ofstream out(savepath);
	nbn_graph.outputNBNnetwork(out);
	out.close();
	ouputNBNdata(saveDri2 + tspname + "_nbn.txt", nbn_graph.get_nbn_data());




}



void filterNBNgap(const std::string& saveDir, const std::string& filename) {
	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;

	inputNBNdata(saveDir + filename + "_nbn.txt", nbn_data);
	ofec::nbn::NBNinfo nbnInfo;
	nbnInfo.fromNBNdata(nbn_data);



}



void runTask_visualization() {
	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	//	ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";

	//ofec::g_working_directory = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";

	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn2/";

	auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn_network2_1/";
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
	readTspInstance(ofec::g_working_directory, filenames);
	//for (auto& it : tasks) {
	//	it.m_tspname = it.m_tspname.substr(0, it.m_tspname.find_first_of("."));
	//}


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

struct GlobalThreadInfo {

	std::list<std::string> m_tspNames;

	std::mutex m_mtx;

	std::string dir1;
	std::string dir2;
	std::string saveDir;
	std::string saveDir2;

};


void runThreadTask(GlobalThreadInfo& info) {

	std::string taskName;

	while (true) {
		taskName.clear();

		{
			info.m_mtx.lock();
			if (!info.m_tspNames.empty()) {
				taskName = info.m_tspNames.front();
				info.m_tspNames.pop_front();

			}
			info.m_mtx.unlock();

		}

		if (!taskName.empty()) {
			auto& it = taskName;
			if (it == "u574.tsp") {
				calculateEAX_LKH_NBNnetwork(info.dir2, info.saveDir, info.saveDir2, it);
			}
			else {
				calculateEAX_LKH_NBNnetwork(info.dir1, info.saveDir, info.saveDir2, it);
			}
		}
	}
}


void runNBNvisualizationMultiTask() {
	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	//	ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";

	ofec::g_working_directory = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";

	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn/";

	auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn_network_v2_local/";
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



	GlobalThreadInfo globalInfo;
	globalInfo.dir1 = dir1;
	globalInfo.dir2 = dir2;
	globalInfo.saveDir = saveDir;
	globalInfo.saveDir2 = saveDir2;


	for (auto& it : filenames) {
		globalInfo.m_tspNames.push_back(it);

	}


	std::vector<std::thread> thrds;
	int num_task = std::thread::hardware_concurrency();
	//	num_task = 1;
	for (size_t i = 0; i < num_task; ++i) {
		thrds.push_back(std::thread(
			&runThreadTask, std::ref(globalInfo)));
	}
	for (auto& thrd : thrds)
		thrd.join();



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


void readTask(const std::string& filepath, std::vector<Info>& tasks) {




	std::ifstream in(filepath);


	in.close();
}


void readTask(std::vector<Info>& tasks) {

	int sizeFrom = tasks.size();
	std::string dir1 = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn_compareResults/tsp_instance2_ver.txt";
	readTxtData(dir1, tasks);
	int sizeTo = tasks.size();
	//for (int idx(sizeFrom); idx < sizeTo; ++idx) {
	//	
	//}

	std::string dir2 = "//172.24.34.11/share/2018/diaoyiya/ofec-data/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn_compareResults/tsp_instance2_ver2.txt";
	readTxtData(dir2, tasks);



	//dir2 = "//172.24.34.11/share/2018/diaoyiya/ofec-data/paper_com_experiment_data/totalTsp/eax_lkh_total_sortedX_nbn_compareResults/tsp_instance2_ver.txt";
	//readTxtData(dir2, tasks);



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

		std::vector<Info> tasks;
		readTask(tasks);

		std::sort(tasks.begin(), tasks.end());

		for (auto& it : tasks) {
			std::cout << it.m_tspname << std::endl;
		}

		int stop = -1;

		//runNBNvisualizationMultiTask();
		//runTask_visualization();
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