
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


// oneMaxInfo

struct OneMaxInfo {
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
}


void analyseOneMax(OneMaxInfo& info) {

	info.nbnInfo.modify();
	

	using namespace ofec::nbn;
	auto& values = info.values;
	auto& nbnInfo = info.nbnInfo;

	values.push_back(calculateNBDstd(nbnInfo));
	values.push_back(calculateMaxNeutralBasins(nbnInfo));

	std::vector<int> optIds;
	getOptimumIds(optIds, nbnInfo);

	values.push_back(optIds.size());

	values.push_back(calculateMaxNBDofOptimum(nbnInfo, optIds));

	std::vector<std::vector<int>> basins;
	calculateBasins(basins, nbnInfo, optIds);

	values.push_back(basinSizesSTD(basins));
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
	std::map< std::pair<std::string, std::pair<int, int>>, OneMaxInfo> totalInfos;
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
		OneMaxInfo oneMaxInfo;
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



void runTask() {
	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	std::string filepath = "yes.txt";
	ofec::inputNBNdata(filepath, nbn_data);

}

namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run(int argc, char* argv[]) {

		using namespace ofec;
		using namespace std;

		registerInstance();
		//runTask();

		std::ofstream out("//172.29.203.176/e/DiaoYiya/paper_com_experiment_data/oneMaxAnalysisData/compareOneMaxData2.txt");

		readFiles("//172.29.203.176/e/DiaoYiya/paper_com_experiment_data/oneMax/onemaxNeighborFilterNetwork2", out);

		out.close();
	}


}