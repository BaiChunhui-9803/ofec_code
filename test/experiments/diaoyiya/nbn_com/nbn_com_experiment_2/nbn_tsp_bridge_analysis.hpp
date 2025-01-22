
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

#include "../utility/visualization/idee_trait.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_info.h"

#include "../utility/matlab/matlab_utility.h"

#include "../utility/nbn_visualization/nbn_fla/tsp_related/nbn_fla_tsp.h"




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



struct NBNanalyseInfo {

	ofec::nbn::NBNinfo nbnInfo;
	std::vector<double> values;


};


void analyseNameOneMax(std::vector<std::string>& values) {

	values.push_back("NBDstd");
	values.push_back("numberOpts");
	values.push_back("MaxNBDofOptimum");
	values.push_back("basinSizesSTD");
	values.push_back("bestNBDstd");


	//values.push_back("lkhNBD");
}


void analyseNBN(NBNanalyseInfo& info, double disRatio) {


	info.nbnInfo.modify();

	using namespace ofec::nbn;
	auto& values = info.values;
	auto& nbnInfo = info.nbnInfo;

	std::vector<int> neurtual_size;

	values.push_back(calculateNBDstd(nbnInfo));
	std::vector<int> optIds;
	getOptimumIdsByDistance(optIds, nbnInfo, disRatio);

	values.push_back(optIds.size());

	values.push_back(calculateMaxNBDofOptimum(nbnInfo, optIds));

	std::vector<std::vector<int>> basins;
	calculateBasins(basins, nbnInfo, optIds);

	values.push_back(basinSizesSTD(basins));
	values.push_back(getNBDvalue(nbnInfo));
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



void readNBNdata(const std::string& saveDir, const std::string& filename, ofec::nbn::NBNinfo& nbnInfo) {
	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	inputNBNdata(saveDir + filename + "_nbn.txt", nbn_data);
	nbnInfo.fromNBNdata(nbn_data);
}



void readNBNinfo(const std::string& saveDir, const std::string& filename, ofec::nbn::NBNinfo& nbnInfo) {
	nbnInfo.inputNBNinfo(saveDir, filename);
	//std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	//inputNBNdata(saveDir + filename + "_nbn.txt", nbn_data);
	//nbnInfo.fromNBNdata(nbn_data);
}




void getNBN_bridgeData(const std::string& filedir, const std::string& saveDir, const std::string& filename) {
	NBNanalyseInfo info;
	readNBNdata(filedir, filename, info.nbnInfo);
	std::vector<int> bridgeIds;
	ofec::nbn::tsp::getNarrowGap(info.nbnInfo, bridgeIds, 0.9, 0.01);


	std::ofstream out(saveDir + filename + "_bridgeIds.txt");
	ofec::matlab::outputVector(out, bridgeIds);
	out.close();
}


void getNBN_bridgeData(const std::string& filedir, const std::string& filename, std::vector<int>& bridgeIds) {
	NBNanalyseInfo info;
	readNBNdata(filedir, filename, info.nbnInfo);
	bridgeIds.clear();
	ofec::nbn::tsp::getNarrowGap(info.nbnInfo, bridgeIds, 0.95, 0.01);


	//std::ofstream out(saveDir + filename + "_bridgeIds.txt");
	//ofec::matlab::outputVector(out, bridgeIds);
	//out.close();
}

void getNBNinfo_bridgeData(const std::string& filedir, const std::string& saveDir, const std::string& filename, NBNanalyseInfo& info, std::vector<int>& bridgeIds) {

	readNBNinfo(filedir, filename, info.nbnInfo);
	//std::vector<int> bridgeIds;
	bridgeIds.clear();


	ofec::nbn::tsp::getNarrowGap(info.nbnInfo, bridgeIds, 0.9, 0.01);



	std::ofstream out(saveDir + filename + "_bridgeIds.txt");
	ofec::matlab::outputVector(out, bridgeIds);
	out.close();
}




void getNBN_bridgeData() {
	const std::string& nbnDataDir = "//172.24.34.11/share/2018/diaoyiya/ofec-data/paper_com_experiment_data/totalTsp/nbn_data_eax_sampling_1000000_final_nbnNetworks_remote/";
	std::string nbnInfoDir = "//172.24.34.11/share/2018/diaoyiya/ofec-data/paper_com_experiment_data/totalTsp/nbn_data_eax_sampling_1000000_final/";
	nbnInfoDir = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread_equalBetter/";
	std::string savedir = "//172.24.34.11/share/2018/diaoyiya/ofec-data/paper_com_experiment_data/totalTsp/nbn_data_eax_sampling_1000000_final_nbnNetworks_bridgeData4/";
	savedir = "//172.24.34.11/share/2018/diaoyiya/ofec-data/paper_com_experiment_data/totalTsp/nbn_sub_region_solVec_hnswModelSingleThread_equalBetter_bridgeData/";


	
	std::filesystem::create_directories(savedir);
	std::vector<std::string> filenames = { "2202","5955", "u574" };
	//	filenames = { "5955"};
	std::vector<int> subRegions = { 3,6,12,25,50,100 };
	//	subRegions = { 6 };
	std::string randomSamples = "_randomSamples";

	std::string neighborK = "_neighborK";

	std::vector<int> bridgeData1;
	//std::vector<int> bridgeData2;



	std::ofstream out(savedir + "BridgeInfo.txt");





	for (auto& curname : filenames) {

		out << curname << "\t";

		NBNanalyseInfo info;
		getNBNinfo_bridgeData(nbnInfoDir, savedir, curname + randomSamples, info, bridgeData1);
		//	getNBN_bridgeData(nbnDataDir, curname + randomSamples, bridgeData2);
		ofec::nbn::tsp::GapInfo gapInfo;
		ofec::nbn::tsp::calculateNarrowGapInfo(info.nbnInfo, bridgeData1, gapInfo);

		out << "origin\t" << gapInfo.m_maxlength << "\t" << gapInfo.m_numNodes << "\t" << gapInfo.m_avgBasinSize << "\t";



		//getNBN_bridgeData(filedir, savedir, curname + randomSamples);
		for (auto& k : subRegions) {



			NBNanalyseInfo info;
			getNBNinfo_bridgeData(nbnInfoDir, savedir, curname + randomSamples + neighborK + "_" + std::to_string(k), info, bridgeData1);
			//	getNBN_bridgeData(nbnDataDir, curname + randomSamples, bridgeData2);
			ofec::nbn::tsp::GapInfo gapInfo;
			ofec::nbn::tsp::calculateNarrowGapInfo(info.nbnInfo, bridgeData1, gapInfo);

			//	getNBNinfo_bridgeData(nbnInfoDir, curname + randomSamples + neighborK + "_" + std::to_string(k), bridgeData1);


			out << k << "\t" << gapInfo.m_maxlength << "\t" << gapInfo.m_numNodes << "\t" << gapInfo.m_avgBasinSize << "\t";

			//getNBN_bridgeData(filedir, savedir, curname + randomSamples + neighborK + "_" + std::to_string(k));
		}
		out << std::endl;
	}

	out.close();
}


namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run(int argc, char* argv[]) {

		using namespace ofec;
		using namespace std;

		registerInstance();

		getNBN_bridgeData();



	}


}