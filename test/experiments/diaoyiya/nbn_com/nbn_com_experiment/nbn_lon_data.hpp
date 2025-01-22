
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
#include "../instance/algorithm/visualize/sampling/instance/sampling_eax_tsp.h"
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

#include "../instance/algorithm/visualize/lon/lon_combinatorial/lon_tsp.h"




// lon related network



void outputLon(std::stringstream& buf, lon::LonInfo& filterLon) {
	UTILITY::vectorToStream(filterLon.m_node_funnel_idxs, buf);
	UTILITY::vectorToStream(filterLon.m_node_funnel_fit, buf);
	UTILITY::vectorToStream(filterLon.m_node_size, buf);
	UTILITY::vectorToStream(filterLon.m_node_fit, buf);
	UTILITY::vvectorToStream(filterLon.m_graph, buf);
}


void inputLon(std::stringstream& buf, lon::LonInfo& lon) {
	UTILITY::streamToVector(lon.m_node_funnel_idxs, buf);
	UTILITY::streamToVector(lon.m_node_funnel_fit, buf);
	UTILITY::streamToVector(lon.m_node_size, buf);
	UTILITY::streamToVector(lon.m_node_fit, buf);
	UTILITY::streamToVvector(lon.m_graph, buf);
	lon.m_nodeId.resize(lon.m_node_fit.size());
	for (int idx(0); idx < lon.m_nodeId.size(); ++idx) {
		lon.m_nodeId[idx] = idx;
	}
}

void outputLon(const std::string& filename, lon::LonInfo& filterLon) {

	std::ofstream out(filename);
	{
		std::stringstream buf;
		outputLon(buf, filterLon);
		out << buf.rdbuf();
	}
	out.close();
}

void inputLon(const std::string& filename, lon::LonInfo& filterLon) {

	std::ifstream in(filename);
	{
		std::stringstream buf;
		buf << in.rdbuf();
		inputLon(buf, filterLon);
	}
	in.close();
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


void inputFromFile(ofec::ParameterVariantStream& paramsStream, const std::string& filepath) {
	std::stringstream buf;
	std::ifstream in(filepath);
	buf << in.rdbuf();
	in.close();
	ofec::variants_stream::stringstream2parameterStream(buf, paramsStream);

}

struct NBNinfo {

	std::vector<std::shared_ptr<ofec::SolutionBase>> m_sols;

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

	static std::string getNBNinfoPath(const std::string& saveDir, const std::string& filename) {
		auto filepath = saveDir + filename + "_nbnInfo.txt";
		return filepath;
	};
	void inputNBNinfo(const std::string& saveDir, const std::string& filename) {
		ofec::ParameterVariantStream paramsStream;
		auto filepath = saveDir + filename + "_nbnInfo.txt";
		inputFromFile(paramsStream, filepath);
		paramsStream >> belong;
		paramsStream >> dis2parent;
		paramsStream >> vFitness;
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
	}
	void inputTspSolutions(const std::string& saveDir, const std::string& filename, ofec::Environment* env) {
		ofec::ParameterVariantStream paramsStream;
		inputFromFile(paramsStream, saveDir + filename + "_solVariables.txt");
		size_t solSize(0);

		paramsStream >> solSize;
		m_sols.resize(solSize);
		auto pro = env->problem();

		for (auto& it : m_sols) {
			it.reset(pro->createSolution());
			auto& cursol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*it);
			paramsStream >> cursol.variable().vect();
		}

		solbases.clear();
		for (auto& it : m_sols) {
			solbases.push_back(it.get());
		}

	}
};



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

void calculateNBNnetwork(const std::string& filepath, NBNinfo& nbnInfo, ofec::Environment* env, ofec::Random* rnd) {

	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(env->problem());

	std::vector<int> solIds(nbnInfo.solbases.size());
	for (int idx(0); idx < solIds.size(); ++idx) {
		solIds[idx] = idx;
	}
	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	ofec::transferNBNData(nbn_data, solIds, nbnInfo.belong, nbnInfo.dis2parent, nbnInfo.vFitness);

	ofec::TreeGraphSimple nbn_graph;
	nbn_graph.setNBNdata(nbn_data);
	nbn_graph.modifyBestSols(rnd, env->problem(), nbnInfo.solbases);
	nbn_graph.calNetwork(tsp_pro->numberVariables());


	std::string savepath = filepath + "_network.txt";
	std::ofstream out(savepath);
	nbn_graph.outputNBNnetwork(out);
	out.close();
	ouputNBNdata(filepath + "_nbn.txt", nbn_graph.get_nbn_data());

}



void calcualteNBNlocalNetwork(NBNinfo& nbnInfo, const std::vector<std::vector<ofec::SolutionBase>>& sols, std::vector<std::pair<int, int>>& solIdAFunIds) {

	struct NBNinfo {
		int m_belong = -1;
		double m_dis2parent = std::numeric_limits<double>::max();
		double m_betterFit = std::numeric_limits<double>::lowest();
	};

	std::vector<std::vector<NBNinfo>> solToNBNId(sols.size());
	for (int idx(0); idx < solToNBNId.size(); ++idx) {
		solToNBNId[idx].resize(sols[idx].size());
		//for (auto& it : solToNBNId[idx]) {
		//	it.first = -1;
		//	it.second = std::numeric_limits<double>::lowest();
		//}
	}

	for (int idx(0); idx < sols.size(); ++idx) {
		for (int idy(0); idy < sols[idx].size(); ++idy) {
			auto& toNBN = solToNBNId[idx][idy];
			auto& cursol = sols[idx][idy];
			for (int idSol(0); idSol < nbnInfo.solbases.size(); ++idSol) {
				auto& otherSol = nbnInfo.solbases[idSol];
				//if(otherSol-)
			}
		}
	}
}


void calTask(const std::string& proReadDir,
	const std::string& lonRealDir,
	const std::string& saveDir,
	const std::string& tspfilename) {
	using namespace ofec;
	using namespace std;
	using namespace chrono;
	auto tspname = tspfilename.substr(0, tspfilename.find_first_of("."));

	std::cout << "calculating filename\t" << tspname << std::endl;

	std::shared_ptr<Environment> env;
	genenrateTSPenv(proReadDir, tspfilename, env);
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

	//std::shared_ptr<ofec::SolutionBase> bestSol;
	lon::LonInfo lon;
	lon::LonInfo filterLon;
	std::vector<std::vector<int>> lonSols;
	std::vector<int> lonSolNodeFit;
	std::vector<int> lonOriginId2newId;


	{

		std::string cursavePath = lonRealDir + tspname;
		cursavePath += "/";


		inputLon(cursavePath + tspname + "_origin_lon.txt", lon);
		lon::filterLON(lon, filterLon, 0.001, lonOriginId2newId);
		ofec::ParameterVariantStream variantStr;

		{
			std::stringstream buf;
			std::ifstream in(cursavePath + tspname + "_sol.txt");
			buf << in.rdbuf();
			in.close();
			ofec::variants_stream::stringstream2parameterStream(buf, variantStr);

		}

		size_t lonSolSize(0);
		variantStr >> lonSolSize;
		lonSols.resize(lonSolSize);
		for (auto& it : lonSols) {
			variantStr >> it;
		}
		//variantStr << lon_sampling.sols();
		variantStr >> lonSolNodeFit;



	}




	NBNinfo nbnInfo;
	for (auto& it : lonSols) {
		if (it.empty())continue;

		std::shared_ptr<ofec::SolutionBase> cursol(pro->createSolution());
		//cursol->initialize(env.get(), rnd);
		auto cursolx = it;
		for (auto& it : cursolx) --it;
		auto& cursolt = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*cursol);
		cursolt.variable().vect() = cursolx;
		//eval_fun(cursolt, env.get());
		nbnInfo.m_sols.push_back(cursol);
	}

	for (auto& it : nbnInfo.m_sols) {
		nbnInfo.solbases.push_back(it.get());
	}

	UTILITY::evaluateRandomSolutionsMultiThreads(nbnInfo.solbases, env.get(), eval_fun);

	//if (nbnInfo.solbases.size() > nbnSamples)
	//	nbnInfo.solbases.resize(nbnSamples);
	std::cout << "calculateNBN solutions by multithread" << std::endl;
	calculateNBN(nbnInfo, env.get(), rnd.get());

	std::cout << "outputNBNinfo solutions by multithread" << std::endl;
	nbnInfo.outputNBNinfo(saveDir, tspname + "_lonSamples");
	nbnInfo.outputTspSolutions(saveDir, tspname + "_lonSamples");


	calculateNBNnetwork(saveDir + tspname + "_lonSamples", nbnInfo, env.get(), rnd.get());


	{
		std::vector<int> nodeFunnelId(lon.m_node_fit.size() - 1, -2);
		for (int idx(1); idx < lon.m_node_fit.size(); ++idx) {
			if (lonOriginId2newId[idx] >= 0) {
				nodeFunnelId[idx - 1] = filterLon.m_node_funnel_idxs[lonOriginId2newId[idx]];
			}
		}


		ofec::ParameterVariantStream paramsStream;
		paramsStream << nodeFunnelId;
		outputToFile(paramsStream, saveDir + tspname + "_lonFunnelIds.txt");
	}
}






void runTask() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
//	ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";


	std::string readDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	readDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_data_eax_sampling_1000000_final/";
	//	saveDir = "F:/DiaoYiya/nbn_com_exp/nbn_sampling/eax_lkh_total_sortedX_2";
	//std::filesystem::create_directories(saveDir);

	readDir = "V:/DiaoYiya/code/data/ofec-data/paper_com_experiment_data/tsp_typical_comparison_exp_read/";


	auto saveDir2 = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/nbn_data_lonSample/";
	std::filesystem::create_directories(saveDir2);

	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames =
	{ "5955.tsp",  "u574.tsp" , "2202.tsp", "1281.tsp" ,"6702.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", };
	filenames =
	{ "5955.tsp",  "u574.tsp" , "2202.tsp", "1281.tsp" };

	//std::reverse(filenames.begin(), filenames.end());

	for (auto& it : filenames) {
		if (it == "u574.tsp") {
			calTask(dir2, readDir, saveDir2, it);
		}
		else {
			calTask(dir1, readDir, saveDir2, it);
		}
	}



	

}



namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run() {

		using namespace ofec;
		using namespace std;

		registerInstance();


		runTask();

	}


}