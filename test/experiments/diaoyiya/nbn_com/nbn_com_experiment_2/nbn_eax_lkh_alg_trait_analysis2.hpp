
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


	values.push_back("lkhNBD");
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
		lkh_alg.run(readDir, tspfilename, totalSols);

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

		it = rnd.uniform.nextNonStd<unsigned>(0, std::numeric_limits<unsigned>::max());
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







void calculateNBN(ofec::nbn::NBNinfo& info, ofec::Environment* env, ofec::Random* rnd) {
	auto& solbases = info.m_solbases;
	nbn::HnswModel hnswModel;
	hnswModel.initialize(env, rnd, std::thread::hardware_concurrency());
	std::vector<int> solIds;
	//	new_datum->m_hnswModel2.setNumberThreads(1);
	hnswModel.addDataMutliThread(solbases, solIds);


	auto& belong = info.m_belong;
	auto& dis2parent = info.m_dis2parent;
	auto& vFitness = info.m_vFitness;
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


void ananlyseNBNtask(
	const std::string& saveDir,
	std::vector<ofec::SolutionBase*>& totalSols,
	NBNanalyseInfo& nbnAnaInfo,
	double disRatio,
	const std::string& taskName,
	ofec::Environment* env, ofec::Random* rnd


) {

	using namespace ofec;
	using namespace std;

	auto pro = env->problem();
	auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
	auto eval_fun =
		[](ofec::SolutionBase& sol, ofec::Environment* env) {
		using namespace ofec;
		sol.evaluate(env, false);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};

	std::vector<SolutionBase*> filterSols;
	std::vector<int> solId2FilterId;
	filterUniqueSols(totalSols, filterSols, solId2FilterId, tsp_pro, rnd);

	UTILITY::evaluateRandomSolutionsMultiThreads(filterSols, env, eval_fun);

	auto& nbnInfo = nbnAnaInfo.nbnInfo;
	nbnInfo.m_solbases = filterSols;
	std::cout << "calculateNBN solutions by multithread" << std::endl;
	calculateNBN(nbnInfo, env, rnd);

	std::cout << "outputNBNinfo solutions by multithread" << std::endl;
	nbnInfo.outputNBNinfo(saveDir, taskName);

	analyseNBN(nbnAnaInfo, disRatio);

}


void runEAX_LKH(const std::string& readDir,
	const std::string& saveDir,
	const std::string& tspfilename,
	NBNanalyseInfo& lkhNBNinfo,
	NBNanalyseInfo& eaxNBNinfo
) {
	using namespace ofec;
	using namespace std;
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
	std::vector<std::vector<std::vector<std::shared_ptr<ofec::SolutionBase>>>> eaxSols;
	runEax(readDir, tspfilename, eaxSols);

	std::vector<std::vector<int>> lkh_ids;
	std::vector<std::vector<std::vector<int>>> eax_ids;

	//std::vector<SolutionBase*> totalSols;
	////	int solId(0);
	//eax_ids.resize(sols.size());
	//for (int idx(0); idx < eax_ids.size(); ++idx) {
	//	eax_ids[idx].resize(sols[idx].size());
	//	for (int idy(0); idy < eax_ids[idx].size(); ++idy) {
	//		eax_ids[idx][idy].resize(sols[idx][idy].size());
	//		for (int idz(0); idz < eax_ids[idx][idy].size(); ++idz) {
	//			eax_ids[idx][idy][idz] = totalSols.size();
	//			totalSols.push_back(sols[idx][idy][idz].get());

	//		}
	//	}
	//}

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
			lkh_ids[idx][idy] = lkh_sols.size();
			lkh_sols.push_back(cursol);
			//	totalSols.push_back(cursol.get());

		}
	}




	std::vector<SolutionBase*> bestSols;
	{
		std::vector<SolutionBase*> totalSols;
		for (auto& it : eaxSols) {
			for (auto& it2 : it.back()) {
				totalSols.push_back(it2.get());
			}
		}
		std::string  name = tspname + "_eax";


		for (auto& it : totalSols) {
			it->setFitness(std::numeric_limits<double>::lowest());
		}

		ananlyseNBNtask(saveDir, totalSols, eaxNBNinfo, 10, name, env.get(), rnd.get());

		double bestFit = std::numeric_limits<double>::lowest();

		for (auto& it : totalSols) {
			if (bestFit < it->fitness()) {
				bestSols.clear();
				bestSols.push_back(it);
			}
			else if (bestFit == it->fitness()) {
				bestSols.push_back(it);
			}
		}
	}

	{
		std::vector<SolutionBase*> totalSols;
		for (auto& it : lkh_sols) {
			totalSols.push_back(it.get());
		}

		for (auto& it : bestSols) {
			totalSols.push_back(it);
		}

		std::string  name = tspname + "_lkh";
		ananlyseNBNtask(saveDir, totalSols, lkhNBNinfo, 0, name, env.get(), rnd.get());
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


void readTspInstance(const std::string& readDir, std::vector<Info>& tasks) {

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



	//tspnames.clear();
	//std::ofstream out(dirpath + "/experiments2.txt");
	for (auto& it : tasks) {
		//out << it.m_tspname << "\t" << it.m_suc_lkh << "\t" << it.m_suc_eax << std::endl;
		std::cout << it.m_tspname << "\t" << it.m_suc_lkh << "\t" << it.m_suc_eax << std::endl;
		//tspnames.push_back(it.m_tspname);
	}
	//out.close();




}



void runTask() {
	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "//172.29.41.69/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "//172.29.204.109/share/Student/2018/YiyaDiao/code_total/data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data";


	std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/eax_lkh_single_nbnInfo/";

	std::filesystem::create_directories(saveDir);

	std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<Info> filenames;
	//{ "2202.tsp","1281.tsp" ,"u574.tsp" ,"6717.tsp" , "6702.tsp" , "7310.tsp", "9225.tsp", "5955.tsp", };
	readTspInstance(ofec::g_working_directory, filenames);




	std::cout << "run EAX vs LKH single run v2" << std::endl;


	//	filenames.clear();
	//	filenames = { "2202.tsp" };

		//	std::reverse(filenames.begin(), filenames.end());


	std::ofstream out(saveDir + "eax_lkh_compare_result.txt");
	{

		std::vector<std::string> headNames11 = { "tspname" ,"suc_lkh" };
		std::vector<std::string> headNames12 = { "suc_eax" };
		std::vector<std::string> headNames2;

		analyseNameOneMax(headNames2);
		for (auto& it : headNames11) {
			out << it << "\t";
		}
		for (auto& it : headNames2) {
			out << it << "\t";
		}



		for (auto& it : headNames12) {
			out << it << "\t";
		}
		for (auto& it : headNames2) {
			out << it << "\t";
		}
		out << std::endl;
	}


	for (auto& it : filenames) {

		NBNanalyseInfo lkhNBNinfo;
		NBNanalyseInfo eaxNBNinfo;
		if (it.m_tspname == "u574.tsp") {
			runEAX_LKH(dir2, saveDir, it.m_tspname, lkhNBNinfo, eaxNBNinfo);
		}
		else {
			runEAX_LKH(dir1, saveDir, it.m_tspname, lkhNBNinfo, eaxNBNinfo);
		}

		out << it.m_tspname << "\t" << it.m_suc_lkh << "\t";

		{
			auto& values = lkhNBNinfo.values;
			for (auto& it2 : values) {
				out << it2 << "\t";
			}
		}

		out << it.m_suc_eax << "\t";

		{
			auto& values = eaxNBNinfo.values;
			for (auto& it2 : values) {
				out << it2 << "\t";
			}
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

		runTask();


	}


}