
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
#include "../instance/algorithm/combination/eax_tsp/eax_tsp_origin/eax_tsp_alg.h"


#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>








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


struct ProblemInfo {


	std::string m_filename;

	std::mutex m_curMtx;

	int m_curTask = 0;
	unsigned m_fromSeed = 1;
	unsigned m_seedStep = 1e3;
	double m_fromSeed2 = 0.01;
	double m_seedStep2 = 0.001;


	//	std::mutex m_leftMtx;
	int m_leftTask = 0;
	std::vector<double> m_lkh_best;
	std::vector<double> m_eax_best;


	void set(const std::string& filename) {
		m_filename = filename;
		m_curTask = 60;
		m_leftTask = 60;
	}


	void getStatisitc(int& numLKHsuc, int& EAXsuc) {
		numLKHsuc = 0;
		EAXsuc = 0;
		double bestValue1(0), bestValue2(0);
		ofec::calMin(m_lkh_best, bestValue1);
		ofec::calMin(m_eax_best, bestValue2);
		bestValue1 = std::min(bestValue1, bestValue2);
		for (auto& it : m_lkh_best) {
			if (it == bestValue1) ++numLKHsuc;
		}
		for (auto& it : m_eax_best) {
			if (it == bestValue1) ++EAXsuc;
		}
	}

};


struct GlobalInfo {

	std::string readDir;

	std::vector<std::unique_ptr<ProblemInfo>> m_totalTask;


	std::ofstream m_out;
	std::mutex m_outMtx;


	std::ofstream m_data_out;


	void output(const std::string& filename, int numLKHsuc, int EAXsuc) {
		m_outMtx.lock();
		m_out << filename << "\t" << numLKHsuc << "\t" << EAXsuc << std::endl;
		std::cout << filename << "\t" << numLKHsuc << "\t" << EAXsuc << std::endl;
		m_outMtx.unlock();

	}


	void output(const std::string& filename, ProblemInfo& info) {
		m_outMtx.lock();
		m_data_out << filename << "\t";
		for (auto& it : info.m_lkh_best) {
			m_data_out << it << "\t";
		}
		for (auto& it : info.m_eax_best) {
			m_data_out << it << "\t";
		}
		m_data_out << std::endl;
		m_outMtx.unlock();

	}
};


double runLKH(const std::string& readDir, const std::string& tspfilename, unsigned seed) {
	std::vector<std::vector<int>> totalSols;
	LKH::LKHAlg lkh_alg;

	lkh_alg.SRandom(seed);
	lkh_alg.run(readDir, tspfilename, totalSols);

	std::shared_ptr<ofec::Environment> env;
	genenrateTSPenv(readDir, tspfilename, env);

	std::vector<int> cursolx;
	lkh_alg.getBestSolution(cursolx);
	for (auto& it : cursolx) --it;
	std::unique_ptr<ofec::SolutionBase> cursol(env->problem()->createSolution());
	auto& tspsol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*cursol);
	tspsol.variable().vect() = cursolx;
	tspsol.evaluate(env.get(), false);
	return tspsol.objective()[0];
}

double runEaxTsp(const std::string& readDir, const std::string& tspfilename, double seed) {
	std::shared_ptr<ofec::Environment> env;
	genenrateTSPenv(readDir, tspfilename, env);
	using namespace ofec;
	env->setAlgorithm(EAX_TSP::create());

	//env->setAlgorithm(Factory<Algorithm>::produce("EAX_TSP"));
//	env->algorithm()->inputParameters().input(localInfo.m_algParam);
	env->algorithm()->recordInputParameters();
	env->initializeAlgorithm(seed);

	env->runAlgorithm();
	auto eax_alg = CAST_EAX_TSP(env->algorithm());
	return eax_alg->bestObjective();

}


void runAlgOneThread(GlobalInfo& globalInfo) {
	int curTask = 0;
	while (curTask < globalInfo.m_totalTask.size()) {
		auto& proInfo = *globalInfo.m_totalTask[curTask];
		bool flagContinue = true;
		while (true) {
			int runId = 0;
			unsigned curSeed1(0);
			double curSeed2(0);
			{

				std::unique_lock lock(proInfo.m_curMtx);
				//proInfo.m_curMtx.lock();
				if (proInfo.m_curTask == 0) {
					flagContinue = false;

					//	break;
				}
				else {
					runId = proInfo.m_curTask--;
					curSeed1 = proInfo.m_fromSeed;
					proInfo.m_fromSeed += proInfo.m_seedStep;
					curSeed2 = proInfo.m_fromSeed2;
					proInfo.m_fromSeed2 += proInfo.m_seedStep2;
					//	proInfo.m_curMtx.unlock();
				}
				//proInfo.m_curMtx.unlock();
				//std::unique_lock lock(proInfo.m_curMtx);


			}
			if (!flagContinue)break;
			double bestValue(0);
			if (runId % 2) {
				//if (runId == 59) {
				//	int stop = -1;
				//}
				bestValue = runLKH(globalInfo.readDir, proInfo.m_filename, curSeed1);
			}
			else {
				bestValue = runEaxTsp(globalInfo.readDir, proInfo.m_filename, curSeed2);
				//	bestValue = 0;
			}

			/*		if (bestValue < 0) {
						int stop = -1;
					}*/
					//bestValue = abs(bestValue);
			{
				proInfo.m_curMtx.lock();
				if (runId % 2) {
					proInfo.m_lkh_best.push_back(bestValue);
				}
				else {
					proInfo.m_eax_best.push_back(bestValue);
				}
				if (--proInfo.m_leftTask == 0) {
					int numLKHsuc, EAXsuc;
					proInfo.getStatisitc(numLKHsuc, EAXsuc);
					globalInfo.output(proInfo.m_filename, numLKHsuc, EAXsuc);
					globalInfo.output(proInfo.m_filename, proInfo);
				}
				proInfo.m_curMtx.unlock();
			}

		}




		++curTask;
	}

}


void readDirFun(const std::string& readDir, std::vector<std::string>& filenames) {

	for (const auto& entry : std::filesystem::directory_iterator(readDir)) {
		//std::cout << entry.path().filename() << std::endl;
		filenames.push_back(entry.path().filename().string());
	}

}


void runTaskMultiThread(const std::string& readDir, const std::string& saveDir) {
	int num_task = std::thread::hardware_concurrency();
	//num_task = 1;
	GlobalInfo globalInfo;
	std::vector<std::string> filenames;
	readDirFun(readDir, filenames);
	//filenames.clear();
	//filenames.push_back("0628.tsp");

	globalInfo.m_totalTask.resize(filenames.size());
	for (int idx(0); idx < globalInfo.m_totalTask.size(); ++idx) {
		globalInfo.m_totalTask[idx].reset(new ProblemInfo);
		globalInfo.m_totalTask[idx]->set(filenames[idx]);
	}


	globalInfo.m_out.open(saveDir + "result.txt");
	globalInfo.m_data_out.open(saveDir + "totalData.txt");
	globalInfo.readDir = readDir;

	std::vector<std::thread> thrds;
	for (size_t i = 0; i < num_task; ++i) {
		thrds.push_back(std::thread(
			runAlgOneThread, std::ref(globalInfo)));
	}
	for (auto& thrd : thrds)
		thrd.join();

	globalInfo.m_data_out.close();
	globalInfo.m_out.close();

}



namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run(int argc, char* argv[]) {


		registerInstance();
		ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
		ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
		ofec::g_working_directory = "//172.29.203.176/e/DiaoYiya/code/data/ofec-data";

		std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
		saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
		saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/lkhEaxComData_test/";

		std::filesystem::create_directories(saveDir);


		std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
		std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";



		runTaskMultiThread(dir1, saveDir);


		//runLKHmultiThread(dir1, saveDir, "2202.tsp");
		//runLKHmultiThread(dir2, saveDir, "u574.tsp");

		//std::vector<std::string> filenames;

		//readDir(dir1, filenames);

		//calTaskMultiThreadInputOutputTest(dir1, saveDir, "2202.tsp");
		//calTaskMultiThreadInputOutputTest(dir2, saveDir, "u574.tsp");




	//	runLKH(dir2, saveDir, "GR666.tsp", 1);


		//runTask(saveDir, dir2, "GR666.tsp");
		//runTask(saveDir, dir1, "2202.tsp");


	//	runLonGen(dir2, saveDir, "GR666.tsp");
	//	runLonGen(dir1, saveDir, "2202.tsp");


		//inputLonInfo(saveDir, "GR666.tsp");

	}


}