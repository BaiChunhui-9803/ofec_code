
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
//#include "../core/algorithm/population.h"
//#include "../utility/nondominated_sorting/filter_sort.h"
//#include "../instance/algorithm/visualize/sampling/instance/sampling_eax_tsp.h"
//#include "../instance/algorithm/visualize/sampling/sampling_multiThread.h"
//#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
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
//#include "../instance/algorithm/combination/eax_tsp/eax_tsp_test/eax_tsp_alg.h"
//#include "../instance/problem/combination/travelling_salesman/tsp_offline_data/travelling_salesman_offline_data.h"

//#include "../utility/nbn_visualization/nbn_fla/nbn_fla.h"
//#include "../instance/algorithm/combination/eax_tsp/eax_tsp_basin/eax_tsp_basin_multipop.h"
//#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
//#include "../instance/algorithm/combination/LKH_origin/INCLUDE/LKH.h"
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include "../core/problem/continuous/continuous.h"







void generateSolutionsConsRadiusSub(
	const ofec::SolutionBase& centerSol,
	double sampleRadius,
	std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, int from, int to, ofec::Environment* env, ofec::Random* rnd) {
	using namespace std;
	using namespace ofec;
	auto conPro = CAST_CONOP(env->problem());
	auto& curCsol = dynamic_cast<const ofec::Continuous::SolutionType&>(centerSol);
	auto boundary = conPro->boundary();

	std::vector<double> radius(boundary.size());
	for (int idx(0); idx < radius.size(); ++idx) {
		radius[idx] = (boundary[idx].second - boundary[idx].first) * sampleRadius / 2.0;;
	}

	auto initboudnary = boundary;
	for (int idx(0); idx < radius.size(); ++idx) {
		initboudnary[idx].first = curCsol.variable()[idx] - radius[idx];
		initboudnary[idx].second = curCsol.variable()[idx] + radius[idx];
	}


	for (int idx(0); idx < radius.size(); ++idx) {
		auto& initRange = initboudnary[idx];
		auto& range = boundary[idx];
		if (initRange.first < range.first) {
			double left = range.first - initRange.first;
			initRange.first = range.first;
			initRange.second += left;
		}
		else if (initRange.second > range.second) {
			double left = initRange.second - range.second;
			initRange.second = range.second;
			initRange.first -= left;
		}
	}


	std::vector<int> data;
	for (int idx(from); idx < to; ++idx) {
		sols[idx].reset(env->problem()->createSolution(centerSol));
		auto& cursol = dynamic_cast<ofec::Continuous::SolutionType&>(*sols[idx]);
		for (int idDim(0); idDim < conPro->numberVariables(); ++idDim) {
			cursol.variable()[idDim] = rnd->uniform.nextNonStd<double>(initboudnary[idDim].first, initboudnary[idDim].second);
		}

	}

}




struct ProInfo {
	// con pro info 
	int m_initPre = 0;
	double m_initRadius = 1.0;

	std::string m_proName;
	std::shared_ptr<ofec::Environment> m_env;
	std::shared_ptr<ofec::Random> m_rnd;
	std::shared_ptr<ofec::SolutionBase> m_center;

	int total_file = 0;
	std::string m_winName = "";
	void initialize(const std::string& winName) {
		m_winName = winName;

		m_initRadius = pow(0.1, m_initPre);

		m_rnd.reset(new ofec::Random(0.5));
		createEnvironment(m_env);
		createCenterSolution();

	}

	std::string getFileName()const {

		return "CON_proName_" + m_proName + "_aroundRadiusPre_" + std::to_string(m_initPre);
	}
	void createEnvironment(std::shared_ptr<ofec::Environment>& env) {
		using namespace std;
		using namespace ofec;
		ofec::ParameterMap params;


		params["problem name"] = m_proName;
		//	params["net-sensor-source case"] = std::string("Net2/sensor_case1/source/case11.txt");
			//	params["dataFile2"] = std::string("case1");
			//	params["dataFile3"] = std::string("case11");
			//	params["use LSTM"] = false;
		//	std::shared_ptr<Environment> env;
		env.reset(Environment::create());
		env->recordInputParameters();
		env->initialize();
		env->setProblem(ofec::Factory<ofec::Problem>::produce(m_proName));
		env->problem()->inputParameters().input(params);
		env->problem()->recordInputParameters();
		env->initializeProblem(0.5);


	}

	void createCenterSolution() {
		using namespace ofec;
		auto pro = m_env->problem();
		auto& optSol = pro->optimaBase()->solutionBase(0);
		m_center.reset(pro->createSolution(optSol));

	}




	void generateSolutionsConRadius(std::vector<std::shared_ptr<ofec::SolutionBase>>& sols) {
		using namespace std;
		using namespace ofec;
		auto env = m_env.get();
		auto rnd = m_rnd.get();
		auto pro = env->problem();
		sols.resize(1e6);

		generateSolutionsConsRadiusSub(
			*m_center, m_initRadius, sols, 0, sols.size(), env, rnd);

	}


};




void evaluateSolsSubTask(std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, int from, int to, ofec::Environment* env, ofec::Random* rnd) {

	auto pro = env->problem();
	ofec::Real pos = (pro->optimizeMode(0) == ofec::OptimizeMode::kMaximize) ? 1 : -1;

	for (int idx(from); idx < to; ++idx) {
		sols[idx]->evaluate(env, false);
		sols[idx]->setFitness(pos * sols[idx]->objective(0));
	}
}

void evaluateSols(std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, ProInfo& info) {
	using namespace std;
	using namespace ofec;

	int num_task = std::thread::hardware_concurrency();
	std::vector<std::shared_ptr<Environment>> envs(num_task);
	std::vector<std::shared_ptr<Random>> rnds(num_task);
	for (auto& it : envs) {
		info.createEnvironment(it);
	}
	auto rnd = info.m_rnd.get();
	for (auto& it : rnds) {
		it.reset(new Random(rnd->uniform.next()));
	}

	std::vector<int> tasks;
	std::vector<std::thread> thrds;

	UTILITY::assignThreads(sols.size(), num_task, tasks);
	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		from_to.first = tasks[i];
		from_to.second = tasks[i + 1];

		thrds.push_back(std::thread(
			evaluateSolsSubTask, std::ref(sols),
			tasks[i], tasks[i + 1], envs[i].get(), rnds[i].get()));
	}
	for (auto& thrd : thrds)
		thrd.join();
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




void calculateNBNsingleThread(NBNinfo& info, ofec::Environment* env, ofec::Random* rnd) {

	auto pro = env->problem();

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


	//for (auto& it : solbases) {
	//	vFitness.push_back(it->fitness());
	//}


	std::vector<std::vector<int>> model_neighbors(solbases.size());
	for (int idx(0); idx < solbases.size(); ++idx) {
		auto& nei = hnswModel.getNeighbors(idx);
		auto& curnei = model_neighbors[idx];
		for (auto& neiId : nei.neighbors()) {
			curnei.push_back(neiId.nodeId());
		}
	}



	//for (int idx(0); idx < solbases.size(); ++idx) {
	//	auto& nei = model_neighbors[idx][0];
	//	auto& curnei = neighbors[idx];
	//	curnei = nei;
	//}

	{
		std::cout << "begin nbn" << std::endl;
		auto start = std::chrono::system_clock::now();
		ofec::NBN_NearestBetterCalculator::calculate(solbases, vFitness, model_neighbors,
			belong, dis2parent,
			pro, rnd);
		auto end = std::chrono::system_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		std::cout << "nbn calculate costs"
			<< double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den
			<< "seconds" << std::endl;
	}
}


void calculateNetwork(NBNinfo& nbnInfo, ofec::TreeGraphSimple& nbn_graph, ofec::Environment* env, ofec::Random* rnd) {
	std::vector<int> solIds(nbnInfo.solbases.size());
	for (int idx(0); idx < solIds.size(); ++idx) {
		solIds[idx] = idx;
	}
	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	ofec::transferNBNData(nbn_data, solIds, nbnInfo.belong, nbnInfo.dis2parent, nbnInfo.vFitness);

	;
	nbn_graph.setNBNdata(nbn_data);
	nbn_graph.modifyBestSols(rnd, env->problem(), nbnInfo.solbases);
	nbn_graph.calNetwork(1.0);
}



void calNBN_iterater(ofec::Environment* pro,
	const std::vector<int>& subArea, int& bestId,
	std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
	std::vector<int>& belong,
	std::vector<double>& dis2parent) {

	if (subArea.empty()) {
		bestId = -1;
		return;
	}

	std::vector<int> sortIdx(subArea);
	std::sort(sortIdx.begin(), sortIdx.end(), [&](int a, int b) {
		if (sols[a]->fitness() == sols[b]->fitness()) {
			return a > b;
		}
		else return sols[a]->fitness() < sols[b]->fitness();
		});
	bestId = sortIdx.back();
	for (int idx(0); idx < sortIdx.size(); ++idx) {
		int ida = sortIdx[idx];
		for (int idy(idx + 1); idy < sortIdx.size(); ++idy) {

			int idb = sortIdx[idy];
			if (sols[ida]->fitness() < sols[idb]->fitness()) {
				double dis = sols[ida]->variableDistance(sols[idb]->variableBase(), pro);
				if (dis < dis2parent[ida]) {
					belong[ida] = idb;
					dis2parent[ida] = dis;
				}
			}
		}
	}

}

void calNBN_iterater_ThreadTask(
	int from, int to,
	ofec::Environment* pro,
	const std::vector<std::vector<int>>& subAreas, std::vector<int>& bestId,
	std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
	std::vector<int>& belong,
	std::vector<double>& dis2parent
) {
	for (int idSub(from); idSub < to; ++idSub) {
		calNBN_iterater(pro, subAreas[idSub], bestId[idSub], sols, belong, dis2parent);
	}
}




void nbnCalculator(
	const std::string& saveDir, const std::string& filename,
	std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, ofec::Environment* env, ofec::Random* rnd) {


	ofec::ParameterVariantStream paramsStream;
	auto filepath = saveDir + filename + "_nbnfit.txt";
	inputFromFile(paramsStream, filepath);
	std::vector<double> vFitness;
	paramsStream >> vFitness;

	NBNinfo nbnInfo;
	nbnInfo.m_sols = sols;
	for (auto& it : nbnInfo.m_sols) {
		nbnInfo.solbases.push_back(it.get());
	}
	nbnInfo.vFitness = vFitness;
	for (int idx(0); idx < nbnInfo.solbases.size(); ++idx) {
		nbnInfo.solbases[idx]->setFitness(vFitness[idx]);
	}
	calculateNBNsingleThread(nbnInfo, env, rnd);
	//calculateNBN(nbnInfo, env, rnd);

	ofec::TreeGraphSimple nbn_graph;
	calculateNetwork(nbnInfo, nbn_graph, env, rnd);

	{
		std::string savepath = saveDir + filename + "_network.txt";
		std::ofstream out(savepath);
		nbn_graph.outputNBNnetwork(out);
		out.close();
		ofec::ouputNBNdata(saveDir + filename + "_nbn.txt", nbn_graph.get_nbn_data());
	}
}


//void nbnCalculatorInt(
//	const std::string& saveDir, const std::string& filename,
//	std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, ofec::Environment* env, ofec::Random* rnd) {
//
//
//	ofec::ParameterVariantStream paramsStream;
//	auto filepath = saveDir + filename + "_nbnfit.txt";
//	inputFromFile(paramsStream, filepath);
//	std::vector<double> vFitness;
//	paramsStream >> vFitness;
//
//	NBNinfo nbnInfo;
//	nbnInfo.m_sols = sols;
//	for (auto& it : nbnInfo.m_sols) {
//		nbnInfo.solbases.push_back(it.get());
//	}
//	nbnInfo.vFitness = vFitness;
//	for (int idx(0); idx < nbnInfo.solbases.size(); ++idx) {
//		nbnInfo.solbases[idx]->setFitness(vFitness[idx]);
//	}
//	calculateNBNsingleThread(nbnInfo, env, rnd);
//	//calculateNBN(nbnInfo, env, rnd);
//
//	ofec::TreeGraphSimple nbn_graph;
//	calculateNetwork(nbnInfo, nbn_graph, env, rnd);
//
//	{
//		std::string savepath = saveDir + filename + "_network.txt";
//		std::ofstream out(savepath);
//		nbn_graph.outputNBNnetwork(out);
//		out.close();
//		ofec::ouputNBNdata(saveDir + filename + "_nbn.txt", nbn_graph.get_nbn_data());
//	}
//}



struct GlobalInfo {
	std::string m_win_name = "winLocal";
	std::list<ProInfo> m_proInfos;
	std::string saveDir;
	std::mutex info_mtx;
};


bool file_exists_with_prefix(std::vector<std::filesystem::path>& results, const std::string& savedir, const std::string& prefix) {
	namespace fs = std::filesystem;
	const std::string extension = ".txt";
	std::error_code ec;
	results.clear();

	// 使用递归目录迭代器遍历当前目录及其子目录
	for (const auto& entry : fs::recursive_directory_iterator(savedir)) {
		if (entry.is_regular_file() && entry.path().filename().string().rfind(prefix, 0) == 0) {
			results.push_back(entry.path().filename());
		}
	}

	// 如果结果不为空，说明找到了至少一个符合条件的文件
	return !results.empty();
}



void runTask(GlobalInfo& info) {

	bool hashTask = true;

	ProInfo localInfo;

	while (hashTask) {
		hashTask = false;
		{
			info.info_mtx.lock();
			if (!info.m_proInfos.empty()) {
				localInfo = info.m_proInfos.front();
				info.m_proInfos.pop_front();
				hashTask = true;
			}
			info.info_mtx.unlock();
		}


		if (hashTask) {
			auto filename = localInfo.getFileName();
			std::vector<std::filesystem::path> results;
			if (file_exists_with_prefix(results, info.saveDir, filename)) {
				continue;
			}
			auto testfilename = filename + info.m_win_name + ".dat";
			{
				std::ofstream out(info.saveDir + testfilename);
				out << "yes" << std::endl;
				out.close();
			}
			// 获取当前时间
			auto now = std::chrono::steady_clock::now();
			// 设置2分钟后的时间
			auto wait_until = now + std::chrono::minutes(2);
			// 让当前线程等待，直到 wait_until
			std::this_thread::sleep_until(wait_until);
			file_exists_with_prefix(results, info.saveDir, filename);
			std::sort(results.begin(), results.end());
			if (results.front() == testfilename) {

				std::cout << "runnig\t" << testfilename << std::endl;
				localInfo.initialize(info.m_win_name);
				std::vector<std::shared_ptr<ofec::SolutionBase>> sols;
				localInfo.generateSolutionsConRadius(sols);

				evaluateSols(sols, localInfo);
				std::vector<double> vFitness;
				for (auto& it : sols) {
					vFitness.push_back(it->fitness());
				}


				ofec::ParameterVariantStream paramsStream;
				paramsStream << vFitness;
				auto filepath = info.saveDir + filename + "_nbnfit.txt";
				outputToFile(paramsStream, filepath);

				nbnCalculator(info.saveDir, filename, sols, localInfo.m_env.get(), localInfo.m_rnd.get());


			}


		}


	}
}



void runTotalTasks() {

	using namespace ofec;
	using namespace std;
	ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";



	ofec::g_working_directory = "//172.29.41.69/share/2018/diaoyiya/ofec-data";
	//ofec::g_working_directory = "//172.29.203.176/e/DiaoYiya/code/data/ofec-data";
	ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data";
	//ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";
	//ofec::g_working_directory = "//172.24.242.8/share/Student/2018/YiyaDiao/code_total/data";


	auto saveDir = ofec::g_working_directory + "/RealworldContinous_total/";

	std::filesystem::create_directories(saveDir);


	GlobalInfo globalInfo;
	globalInfo.saveDir = saveDir;
	globalInfo.m_win_name = "_winServer_1";
	globalInfo.m_win_name = "winLocal";
	//std::list<ProInfo>infos;
	auto& infos = globalInfo.m_proInfos;
	std::vector<std::string> proNames = { "gear-train", "FM-sound-waves", "GTOP_Cassini2" ,"RW_CEC2011_T2" };

	std::vector<int> precisions = { 0,1,2,3,4,5,6,7 };


	//double curConRatio = 0.01;
	for (auto& proName : proNames) {
		for (auto& curphase : precisions) {
			ProInfo localInfo;
			localInfo.m_initPre = curphase;
			localInfo.m_proName = proName;
			infos.push_back(localInfo);
		}
	}


	runTask(globalInfo);


}

namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run() {

		using namespace ofec;
		using namespace std;

		registerInstance();
		//runTask();
		runTotalTasks();


	}


}