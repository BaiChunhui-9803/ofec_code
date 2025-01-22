
#include "../core/global.h"
#include "../core/exception.h"
#include <fstream>
#include <list>
#include <iostream>
#include <chrono>
#include <thread>
//#include <nbnFigureOut.h>
#include <iostream>
//#include "testOFECmatlab.h"
//#include "drawMeshFun.h"



#ifdef  OFEC_PYTHON
#include "../utility/python/python_caller.h"
#include "Python.h"
#endif //  OFEC_PYTHON



#include "../utility/nbn_visualization/tree_graph_simple.h"
#include "../utility/nbn_visualization/nbn_onemax_simple.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla.h"

#include "../instance/problem/combination/wmodel/wmodel.h"
#include "../core/problem/continuous/continuous.h"


#include "../utility/general_multithread/general_multithread.h"
#include "../utility/ioh_experimenter/include/ioh/problem/wmodel/wmodel_one_max.hpp"
#include "../utility/ioh_experimenter/include/ioh/problem/wmodel/wmodel_leading_ones.hpp"





void judgeNBN(const std::string& inputDir,
	const std::string& outputDir,
	const std::string filename//, 
/*	std::ofstream& err_out, std::mutex & out_mutx*/) {
	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	std::string filepath = inputDir + filename;
	ofec::inputNBNdata(filepath, nbn_data);
	int numError(0);
	std::vector<int> error_idx;
	int numbelong1 = 0;
	for (auto& it : nbn_data) {
		if (it.m_belong != -1) {
			if (it.m_fitness > nbn_data[it.m_belong].m_fitness) {
				++numError;
				error_idx.push_back(it.m_id);
			}
		}
		else ++numbelong1;
	}
	if (numError == 1 && numbelong1 == 0) {
		std::cout << "ok\t" << filepath << std::endl;
		nbn_data[error_idx.front()].m_belong = nbn_data[error_idx.front()].m_id;
		nbn_data[error_idx.front()].m_dis2parent = 1.0;

		ofec::ouputNBNdata(outputDir + filename, nbn_data);
	}
	else if (numError == 0) {
		std::cout << "ok\t" << filepath << std::endl;

		ofec::ouputNBNdata(outputDir + filename, nbn_data);
	}
	else {
		std::string error_info = "error at \t" + filename + "numError\t" + std::to_string(numError);
		throw ofec::Exception(error_info);
		//std::lock_guard<std::mutex> guard(out_mutx);
		//err_out << filepath << std::endl;
	}
}


class ThreadInfoNBN :public ofec::GeneralMultiThreadInfo {
public:
	//std::ofstream& m_errOut;
	//std::mutex& m_errOut_mtx;
	std::string m_filename;
	std::string m_readfilename;
	size_t m_dim = 0;
	std::string m_proname;
	int m_curRunId = 0;
	std::string m_save_dir;
	bool m_flagOpt = false;


	int m_numSample = 1e6;



	// for wmodel

	int m_instance_id = 0;
	double m_dummy_select_rate = 0.0;
	int m_epistasis_block_size = 0;
	int m_neutrality_mu = 0;
	int m_ruggedness_gamma = 0;
	double m_rug_ratio = 0;


	std::string getConTaskName() {
		return m_proname + "_Dim_" + std::to_string(m_dim);
	}


	std::string getProblemTotalName() {
		std::string proFileName = "OneMax_" "_Dim_" + std::to_string(m_dim)
			+ "_InstanceId_" + std::to_string(m_instance_id)
			+ "_DummySelectRate_" + std::to_string(int(m_dummy_select_rate * 1000))
			+ "_EpistasisBlockSize_" + std::to_string(m_epistasis_block_size)
			+ "_NeutralityMu_" + std::to_string(m_neutrality_mu)
			+ "_RuggednessGamma_" + std::to_string(m_ruggedness_gamma)
			;
		return std::move(proFileName);
	}
};
class ThreadGlobalNBN :public ofec::GeneralMultiThreadInfo {
public:
	std::string m_read_dir;
	std::string m_save_dir;
	std::string m_nbn_subfix = "_nbn.txt";
	int m_totalRun = 0;
	//	std::string m_save_path;
	int m_maxSample = 2e6;


	std::ofstream m_out;
	std::mutex m_mtx;


	void outputHead() {
		m_out << "instance_id\tdummy_select_rate\tepistasis_block_size\t";
		m_out << "neutrality_mu\t<<m_ruggedness_gamma<<\t";
	}

	void output(ThreadInfoNBN& cur) {
		m_out << cur.m_instance_id << "\t" << cur.m_dummy_select_rate << "\t" << cur.m_epistasis_block_size << "\t";
		m_out << cur.m_neutrality_mu << "\t" << cur.m_rug_ratio << "\t";
	}
};


//std::vector<task> m_task;



void testFun()
{
	const std::vector<int> x1(16, 1);
	const std::vector<int> x0(16, 0);

	auto problem_om1 = ioh::problem::wmodel::WModelOneMax(1, 16, 0, 0, 0, 0);
	double value = problem_om1(x1);


	int stop = -1;
	//EXPECT_DOUBLE_EQ(problem_om1(x1), 16.0);
	//EXPECT_DOUBLE_EQ(problem_om1(x0), 0);

	auto problem_lo1 = ioh::problem::wmodel::WModelLeadingOnes(1, 16, 0, 0, 0, 0);
	//EXPECT_DOUBLE_EQ(problem_lo1(x1), 16.0);
	//EXPECT_DOUBLE_EQ(problem_lo1(x0), 0);

	auto problem_om2 = ioh::problem::wmodel::WModelOneMax(1, 16, 0.5, 0, 0, 0);
	//EXPECT_DOUBLE_EQ(problem_om2(x1), 8);
	//EXPECT_DOUBLE_EQ(problem_om2(x0), 0);

	auto problem_lo2 = ioh::problem::wmodel::WModelLeadingOnes(1, 16, 0.5, 0, 0, 0);
	//	EXPECT_DOUBLE_EQ(problem_lo2(x1), 8);
	///	EXPECT_DOUBLE_EQ(problem_lo2(x0), 0);

	auto problem_om3 = ioh::problem::wmodel::WModelOneMax(1, 16, 0., 0, 0, 0);
	//	EXPECT_DOUBLE_EQ(problem_om3(x1), 16);
	//	EXPECT_DOUBLE_EQ(problem_om3(x0), 0);

	auto problem_lo3 = ioh::problem::wmodel::WModelLeadingOnes(1, 16, 0, 0, 0, 0);
	//	EXPECT_DOUBLE_EQ(problem_lo3(x1), 16);
	//	EXPECT_DOUBLE_EQ(problem_lo3(x0), 0);
}





void runOneMax(ThreadGlobalNBN& curGlobalInfo,
	ThreadInfoNBN& curcurInfo) {
	//auto& curGlobalInfo = dynamic_cast<const ThreadGlobalNBN&>(*globalInfo);
	//auto& curcurInfo = dynamic_cast<ThreadInfoNBN&>(*curInfo);
	std::string dir = curGlobalInfo.m_save_dir;

	using namespace ofec;
	ParameterMap params;
	params["problem name"] = std::string("OneMaxWModel");
	params["number of variables"] = curcurInfo.m_dim;



	std::shared_ptr<Environment> env(Environment::create());
	env->recordInputParameters();
	env->initialize();

	env->setProblem(OneMaxWModel::create());
	env->problem()->inputParameters().input(params);
	//env->problem()->inputParameters().at("number of variables")->setValue(curcurInfo.m_dim);
	env->problem()->recordInputParameters();


	auto pro = env->problem();
	auto onemaxPro = GET_OneMaxWModel(pro);


	auto rnd = std::shared_ptr<Random>(new Random(0.5));;


	GET_OneMaxWModel(pro)->setInstanceId(curcurInfo.m_instance_id);
	GET_OneMaxWModel(pro)->setNeutralityMu(curcurInfo.m_neutrality_mu);
	GET_OneMaxWModel(pro)->setDummySelectRate(curcurInfo.m_dummy_select_rate);
	GET_OneMaxWModel(pro)->setEpistasisBlockSize(curcurInfo.m_epistasis_block_size);
	GET_OneMaxWModel(pro)->setRuggednessGamma(curcurInfo.m_ruggedness_gamma);



	pro->initialize(env.get());

	auto eval_fun =
		[](SolutionBase& sol, Environment* env) {
		using namespace ofec;
		auto pro = env->problem();
		sol.evaluate(env, false);
		ofec::Real pos = (pro->optimizeMode(0) == ofec::OptimizeMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};


	curcurInfo.m_error_name = GET_OneMaxWModel(pro)->getProblemTotalName();


	int numSamples = curcurInfo.m_numSample;
	// for test
//	numSamples = 1000;

	std::vector<std::shared_ptr<SolutionBase>> sols;
	std::vector<int> solIds;

	if (onemaxPro->numberVariables() < 20) {
		unsigned long long value = 0;
		unsigned long long totalNum = pow(2, onemaxPro->numberVariables());
		//	std::vector<int> solVec(pro->numberVariables(),0);		
		std::shared_ptr<SolutionBase> curSol;
		for (unsigned long long idx(0); idx < totalNum; ++idx) {
			curSol.reset(pro->createSolution());
			auto& curSolOneMax(dynamic_cast<ProblemVariableVector<int>::SolutionType&>(*curSol));
			unsigned long long bitOne = 1;
			for (int idDim(0); idDim < onemaxPro->numberVariables(); ++idDim) {
				if (idx & bitOne) {
					curSolOneMax.variable()[idx] = 1;
				}
				else curSolOneMax.variable()[idx] = 0;
				bitOne = bitOne << 1;
			}

			curSolOneMax.evaluate(env.get());
			sols.push_back(curSol);
		}

	}
	else {



		{
			std::vector<bool> solution_tag(pro->numberVariables());
			std::set<std::vector<bool>> total_sols;
			std::shared_ptr<SolutionBase> curSol;

			int maxLoop = 5e6;
			while (true) {
				curSol.reset(pro->createSolution());
				auto& curSolOneMax(dynamic_cast<OneMax::SolutionType&>(*curSol));
				while (true) {
					curSolOneMax.initialize(env.get(), rnd.get());
					std::fill(solution_tag.begin(), solution_tag.end(), false);
					for (int idx(0); idx < solution_tag.size(); ++idx) {
						if (curSolOneMax.variable()[idx] == 1) {
							solution_tag[idx] = true;
						}
					}
					--maxLoop;
					if (total_sols.find(solution_tag) == total_sols.end()) {
						total_sols.insert(solution_tag);
						break;
					}
					if (maxLoop == 0)break;
				}
				if (maxLoop == 0)break;

				curSolOneMax.evaluate(env.get());
				sols.push_back(curSol);
				if (sols.size() == numSamples)break;

			}
		}
		//{

		//	const std::string dir = "//172.24.24.151/e/DiaoYiya/experiment_data/onemax/";

		//	std::string filename = GET_OneMaxWModel(pro)->getProblemTotalName() + "_numSamples_" + std::to_string(solIds.size());
		//	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
		//	ofec::inputNBNdata(dir + filename+"_nbn.txt", nbn_data);
		//	ofec::TreeGraphSimple nbn_graph;
		//	nbn_graph.m_task_name = GET_OneMaxWModel(pro)->getProblemTotalName();

		//	nbn_graph.setNBNdata(nbn_data);
		//	nbn_graph.calNetwork();
		//	std::ofstream out(dir + filename + "_network.txt");
		//	nbn_graph.outputNBNnetwork(out);
		//	out.close();
		//}
	}



	solIds.resize(sols.size());
	for (int idx(0); idx < solIds.size(); ++idx) {
		solIds[idx] = idx;
	}


	ofec::NBN_OneMaxSimpleDivison solDivision;
	solDivision.initialize(env, rnd, eval_fun);
	solDivision.setSol(sols, solIds);
	int numThread = std::thread::hardware_concurrency();
	solDivision.setNumTaskUpdateNetwork(numThread);
	solDivision.setNumberThread(numThread);
	solDivision.setNumberLoop(30);

	std::vector<double> belongDis;
	std::vector<int> belong;
	solDivision.initNBN();

	solDivision.udpateNBNmultithread();
	solDivision.getResult(belong, belongDis);
	std::vector<double> fitness;
	solDivision.getFitness(fitness);

	std::string filename = GET_OneMaxWModel(pro)->getProblemTotalName() + "_numSamples_" + std::to_string(solIds.size());





	//// calculate FLA
	//{
	//	ofec::NBN_FLA_info info;
	//	info.m_funId = 1;
	//	info.m_instanceId = 1;
	//	info.m_ruggness_eps = 0.01;
	//	info.calculate(pro, rnd.get(), sols, fitness, belong, belongDis);


	//	double rug = info.m_nbd_rug;
	//	double neutrality = info.m_maxNeutural;
	//	{
	//		std::unique_lock<std::mutex> lock(curGlobalInfo.m_mtx);
	//		curGlobalInfo.output(curcurInfo);
	//		curGlobalInfo.m_out << rug << "\t" << neutrality << "\t";
	//		curGlobalInfo.m_out << std::endl;

	//		//std::cout << filename << "\t" << "is finished" << std::endl;
	//	}

	//}



	{
		std::ofstream nbnOut(dir + filename + "_nbn.txt");// (globalInfo.m_save_dir + filepath + "_nbn.txt");
		ofec::outputNBNdata(nbnOut, solIds, belong, belongDis, fitness);
		nbnOut.close();

	}
	//{

	//	std::vector<ofec::TreeGraphSimple::NBNdata> nbn_data;
	//	nbn_data.resize(solIds.size());
	//	for (int idx(0); idx < nbn_data.size(); ++idx) {
	//		nbn_data[idx].m_id = nbn_data[idx].m_originId = solIds[idx];
	//		nbn_data[idx].m_belong = belong[idx];
	//		nbn_data[idx].m_dis2parent = belongDis[idx];
	//		nbn_data[idx].m_fitness = fitness[idx];

	//	}
	//	//ofec::filterNBNDataByDensity(outputDir, filename, nbn_data);
	//	ofec::TreeGraphSimple nbn_graph;
	//	nbn_graph.m_task_name = GET_OneMaxWModel(pro)->getProblemTotalName();

	//	nbn_graph.setNBNdata(nbn_data);
	//	nbn_graph.calNetwork();
	//	std::ofstream out(dir + filename + "_network.txt");
	//	nbn_graph.outputNBNnetwork(out);
	//	out.close();
	//}




}

void runOneMaxTotalTask() {
	std::vector<int> dims = { 40,80,120,160 ,200 };
	std::vector<int> ruggedness[2] = { {0,1, 2,4,6,8 ,10 } ,{0,1,4} };
	std::vector<int> epistasis[2] = { {0,1,2,3,4,5,6,7,8,9,10,20,31} ,{0,1,2,3,4} };
	std::vector<int> neutral[2] = { { 0,1,2,3,5,7,9,10,20,30 },{0,1,4} };

	//	std::vector<std::vector<int>>  totalIds = { {1,0,0}, };



	ThreadGlobalNBN curGlobalInfo;
	ThreadInfoNBN curcurInfo;


	curGlobalInfo.m_save_dir = "//172.24.24.151/e/DiaoYiya/experiment_data/onemaxtest/";


	curGlobalInfo.m_save_dir = "/mnt/share151e/DiaoYiya/experiment_data/onemax_total/";
	std::filesystem::create_directories(curGlobalInfo.m_save_dir);


	std::ofstream error_out(curGlobalInfo.m_save_dir + "_error_info.txt");



	curcurInfo.m_instance_id = 1;

	int taskId = 0;

	std::vector<int> cur_idxs = { 0,0,0 };
	for (int idT(0); idT < cur_idxs.size(); ++idT) {
		cur_idxs[idT] = 1;

		for (int idD = 0; idD < dims.size(); ++idD) {
			curcurInfo.m_dim = dims[idD];
			for (int idR = 0; idR < ruggedness[cur_idxs[0]].size(); ++idR) {
				curcurInfo.m_ruggedness_gamma = ruggedness[cur_idxs[0]][idR];
				for (int idE(0); idE < epistasis[cur_idxs[1]].size(); ++idE) {
					curcurInfo.m_epistasis_block_size = epistasis[cur_idxs[1]][idE];
					for (int idN(0); idN < neutral[cur_idxs[2]].size(); ++idN) {
						curcurInfo.m_neutrality_mu = neutral[cur_idxs[2]][idN];

						std::cout << "running taskId\t" << taskId << std::endl;
						try {
							clock_t start, end; //定义clock_t变量
							start = clock(); //开始时间
							runOneMax(curGlobalInfo, curcurInfo);

							end = clock();   //结束时间
							std::cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;  //输出时间（单位：ｓ）

						}
						catch (const std::exception& e) {
							std::string info = std::string(e.what()) + " at " +
								curcurInfo.m_error_name;
							error_out << info << std::endl;
							std::cout << "error_info\t" << info << std::endl;
						}
						catch (const ofec::Exception& e) {
							std::string info = std::string(e.what()) + " at " +
								curcurInfo.m_error_name;
							error_out << info << std::endl;
							std::cout << "error_info\t" << info << std::endl;
						}




						std::cout << "finish taskId\t" << curcurInfo.m_error_name << std::endl;



					}
				}
			}
		}


		cur_idxs[idT] = 0;
	}


	error_out.close();


}


void calTask(std::unique_ptr<ofec::GeneralMultiThreadInfo>& curInfo,
	std::unique_ptr<ofec::GeneralMultiThreadInfo>& globalInfo) {
	auto& curGlobalInfo = dynamic_cast<ThreadGlobalNBN&>(*globalInfo);
	auto& curcurInfo = dynamic_cast<ThreadInfoNBN&>(*curInfo);

	runOneMax(curGlobalInfo, curcurInfo);
}


void setTask(ofec::GeneralMultiThread& curThread) {

	ThreadInfoNBN curcurInfo;
	int sampleSize = 1e6;
	//sampleSize = 10000;

	curcurInfo.m_numSample = sampleSize;
	curcurInfo.m_instance_id = 1;
	curcurInfo.m_dim = 80;
	curcurInfo.m_ruggedness_gamma = 0;
	curcurInfo.m_epistasis_block_size = 0;
	curcurInfo.m_neutrality_mu = 0;



	int taskId = 0;
	for (int rugPar = 0; rugPar <= 20; ++rugPar) {
		curcurInfo.m_ruggedness_gamma = rugPar;
		curcurInfo.m_error_name = curcurInfo.getProblemTotalName();
		curThread.ms_info_buffer.emplace_back(new ThreadInfoNBN(curcurInfo));
	}
	curcurInfo.m_ruggedness_gamma = 0;

	for (int rugPar = 0; rugPar <= 30; ++rugPar) {
		curcurInfo.m_neutrality_mu = rugPar;
		curcurInfo.m_error_name = curcurInfo.getProblemTotalName();
		curThread.ms_info_buffer.emplace_back(new ThreadInfoNBN(curcurInfo));
	}



}


void runMultiTask() {
	ofec::GeneralMultiThread curThread;
	curThread.ms_cout_flag = true;
	curThread.ms_global_info.reset(new ThreadGlobalNBN);
	ThreadGlobalNBN& curGlobalInfo = dynamic_cast<ThreadGlobalNBN&>(*curThread.ms_global_info);



	curGlobalInfo.m_save_dir = "//172.24.24.151/e/DiaoYiya/experiment_data/onemaxtest/";
	curGlobalInfo.m_save_dir = "/mnt/share151/DiaoYiya/experiment_data/onemax_nbn_fla/";
	curGlobalInfo.m_save_dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/nbn_fla_onemax/";
	std::filesystem::create_directories(curGlobalInfo.m_save_dir);
	std::cout << "create dir\t" << curGlobalInfo.m_save_dir << std::endl;
	std::ofstream error_out(curGlobalInfo.m_save_dir + "_error_info.txt");



	//		totalInfo.m_save_dir = "//172.24.24.151/e/DiaoYiya/experiment_data/mmop_multiOpts_lon_left/";
	bool m_is_multiTask = true;
	curThread.setNumberThread(80);

	curThread.ms_fun = std::bind(calTask, std::placeholders::_1, std::placeholders::_2);
	setTask(curThread);


	curGlobalInfo.m_out.open(curGlobalInfo.m_save_dir + "finalResult.txt");
	curGlobalInfo.outputHead();
	curGlobalInfo.m_out << "ruggness" << "\t" << "neutrality" << "\t";
	curGlobalInfo.m_out << std::endl;

	std::ofstream out(curGlobalInfo.m_save_dir + "error.txt");
	curThread.run(out);
	out.close();
}


void runOneMaxTotalTaskAnalysis() {
	//std::vector<int> dims = { 80};
	//std::vector<int> ruggedness[2] = { {0,1, 2,4,6,8 ,10 } ,{0,1,4} };
	//std::vector<int> epistasis=  {0,1,2,3,4,5,6,7,8,9,10,20,31} ,{0,1,2,3,4} };
	//std::vector<int> neutral[2] = { { 0,1,2,3,5,7,9,10,20,30 },{0,1,4} };

	//	std::vector<std::vector<int>>  totalIds = { {1,0,0}, };



	ThreadGlobalNBN curGlobalInfo;
	ThreadInfoNBN curcurInfo;


	curGlobalInfo.m_save_dir = "//172.24.24.151/e/DiaoYiya/experiment_data/onemaxtest/";
	curGlobalInfo.m_save_dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data/nbn_visualization_onemax/";
	std::filesystem::create_directories(curGlobalInfo.m_save_dir);
	std::ofstream error_out(curGlobalInfo.m_save_dir + "_error_info.txt");

	curGlobalInfo.m_out.open(curGlobalInfo.m_save_dir + "finalResult.txt");
	curGlobalInfo.outputHead();
	curGlobalInfo.m_out << "ruggness" << "\t" << "neutrality" << "\t";
	curGlobalInfo.m_out << std::endl;


	int max_neutraility_mu = 60;
	int maxRugLoop = 5;
	int maxNeuLoop = 5;
	int maxELoop = 5;
	std::vector<int> dims = { 120,160,40,80 };

	curcurInfo.m_instance_id = 1;
	curcurInfo.m_dim = 80;
	curcurInfo.m_ruggedness_gamma = 0;
	curcurInfo.m_epistasis_block_size = 0;
	curcurInfo.m_neutrality_mu = 0;
	curcurInfo.m_rug_ratio = 0;





	int taskId = 0;


	for (int curdim : dims) {

		curcurInfo.m_dim = curdim;
		curcurInfo.m_ruggedness_gamma = 0;
		curcurInfo.m_rug_ratio = 0;
		curcurInfo.m_neutrality_mu = 0;
		curcurInfo.m_epistasis_block_size = 0;

		for (int rugPar = 0; rugPar <= maxRugLoop; ++rugPar) {

			double ratio = double(rugPar) / double(maxRugLoop);
			curcurInfo.m_rug_ratio = ratio;
			curcurInfo.m_ruggedness_gamma = curcurInfo.m_dim * (curcurInfo.m_dim + 1) / 2.0 * ratio;
			std::cout << "running taskId\t" << taskId++ << std::endl;
			try {
				clock_t start, end; //定义clock_t变量
				start = clock(); //开始时间
				runOneMax(curGlobalInfo, curcurInfo);

				end = clock();   //结束时间
				std::cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;  //输出时间（单位：ｓ）
				std::cout << "curTask\t" << curcurInfo.getProblemTotalName() << std::endl;
			}
			catch (const std::exception& e) {
				std::string info = std::string(e.what()) + " at " +
					curcurInfo.m_error_name;
				error_out << info << std::endl;
				std::cout << "error_info\t" << info << std::endl;
			}
			catch (const ofec::Exception& e) {
				std::string info = std::string(e.what()) + " at " +
					curcurInfo.m_error_name;
				error_out << info << std::endl;
				std::cout << "error_info\t" << info << std::endl;
			}
		}


		curcurInfo.m_ruggedness_gamma = 0;
		curcurInfo.m_rug_ratio = 0;
		curcurInfo.m_neutrality_mu = 0;
		curcurInfo.m_epistasis_block_size = 0;

		for (int rugPar = 0; rugPar <= maxNeuLoop; ++rugPar) {

			double ratio = double(rugPar) / double(maxNeuLoop);

			curcurInfo.m_neutrality_mu = rugPar * max_neutraility_mu / maxNeuLoop;
			std::cout << "running taskId\t" << taskId++ << std::endl;
			try {
				clock_t start, end; //定义clock_t变量
				start = clock(); //开始时间
				runOneMax(curGlobalInfo, curcurInfo);

				end = clock();   //结束时间
				std::cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;  //输出时间（单位：ｓ）
				std::cout << "curTask\t" << curcurInfo.getProblemTotalName() << std::endl;
			}
			catch (const std::exception& e) {
				std::string info = std::string(e.what()) + " at " +
					curcurInfo.m_error_name;
				error_out << info << std::endl;
				std::cout << "error_info\t" << info << std::endl;
			}
			catch (const ofec::Exception& e) {
				std::string info = std::string(e.what()) + " at " +
					curcurInfo.m_error_name;
				error_out << info << std::endl;
				std::cout << "error_info\t" << info << std::endl;
			}
		}



		curcurInfo.m_ruggedness_gamma = 0;
		curcurInfo.m_rug_ratio = 0;
		curcurInfo.m_neutrality_mu = 0;
		curcurInfo.m_epistasis_block_size = 0;

		for (int rugPar = 0; rugPar <= maxELoop; ++rugPar) {


			curcurInfo.m_epistasis_block_size = rugPar * 4 + 2;
			std::cout << "running taskId\t" << taskId++ << std::endl;
			try {
				clock_t start, end; //定义clock_t变量
				start = clock(); //开始时间
				runOneMax(curGlobalInfo, curcurInfo);

				end = clock();   //结束时间
				std::cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;  //输出时间（单位：ｓ）
				std::cout << "curTask\t" << curcurInfo.getProblemTotalName() << std::endl;
			}
			catch (const std::exception& e) {
				std::string info = std::string(e.what()) + " at " +
					curcurInfo.m_error_name;
				error_out << info << std::endl;
				std::cout << "error_info\t" << info << std::endl;
			}
			catch (const ofec::Exception& e) {
				std::string info = std::string(e.what()) + " at " +
					curcurInfo.m_error_name;
				error_out << info << std::endl;
				std::cout << "error_info\t" << info << std::endl;
			}
		}

	}


	error_out.close();


}




void generateSolsRandom(ofec::Environment* env, ofec::Random* rnd,
	std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, int numSols) {
	using namespace ofec;
	sols.resize(numSols * 2);

	auto pro = env->problem();

	for (auto& it : sols) {
		it.reset(pro->createSolution());
		it->initialize(env, rnd);
	}
	std::vector<bool> flag(sols.size(), false);

	std::vector<const ofec::SolutionBase*> sol_ptrs(sols.size());
	for (int idx(0); idx < sols.size(); ++idx) {
		sol_ptrs[idx] = sols[idx].get();
	}
	CAST_ONEMAX(pro)->filterSameSols(sol_ptrs, flag);

	std::vector<std::shared_ptr<ofec::SolutionBase>> filterSols;
	for (int id(0); id < sols.size(); ++id) {
		if (flag[id]) {
			filterSols.push_back(sols[id]);
		}
	}
	if (filterSols.size() > numSols) {
		filterSols.resize(numSols);
	}
	swap(sols, filterSols);


	for (auto& it : sols) {
		it->evaluate(env, false);
	}
}


void generateSols(ofec::Environment* env, ofec::Random* rnd, ofec::SolutionBase& centerSol,
	std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, int neighborK, int numSols) {
	using namespace ofec;
	sols.resize(numSols * 2);

	auto pro = env->problem();

	for (auto& it : sols) {
		it.reset(pro->createSolution());
		CAST_ONEMAX(pro)->initSolutionNearBy(centerSol, *it, neighborK, rnd);
	}

	std::vector<bool> flag(sols.size(), false);

	std::vector<const ofec::SolutionBase*> sol_ptrs(sols.size());
	for (int idx(0); idx < sols.size(); ++idx) {
		sol_ptrs[idx] = sols[idx].get();
	}
	CAST_ONEMAX(pro)->filterSameSols(sol_ptrs, flag);

	std::vector<std::shared_ptr<ofec::SolutionBase>> filterSols;
	for (int id(0); id < sols.size(); ++id) {
		if (flag[id]) {
			filterSols.push_back(sols[id]);
		}
	}
	if (filterSols.size() > numSols) {
		filterSols.resize(numSols);
	}
	swap(sols, filterSols);


	for (auto& it : sols) {
		it->evaluate(env, false);
	}
}



void calNBNonexMaxNeighbor(ofec::Environment* env, ofec::Random* rnd,
	const std::string& dir,
	int neighborK,
	std::vector<std::shared_ptr<ofec::SolutionBase>>& sols) {

	using namespace ofec;
	auto pro = env->problem();
	std::vector<int> solIds;
	solIds.resize(sols.size());
	for (int idx(0); idx < solIds.size(); ++idx) {
		solIds[idx] = idx;
	}

	auto eval_fun =
		[](SolutionBase& sol, Environment* env) {
		using namespace ofec;
		auto pro = env->problem();
		sol.evaluate(env, false);
		ofec::Real pos = (pro->optimizeMode(0) == ofec::OptimizeMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};


	ofec::NBN_OneMaxSimpleDivison solDivision;
	solDivision.initialize(env->getSharedPtr(), rnd->getSharedPtr(), eval_fun);
	solDivision.setSol(sols, solIds);
	int numThread = std::thread::hardware_concurrency();
	solDivision.setNumTaskUpdateNetwork(numThread);
	solDivision.setNumberThread(numThread);
	solDivision.setNumberLoop(30);

	std::vector<double> belongDis;
	std::vector<int> belong;
	solDivision.initNBN();

	solDivision.udpateNBNmultithread();
	solDivision.getResult(belong, belongDis);
	std::vector<double> fitness;
	solDivision.getFitness(fitness);

	std::string filename = GET_OneMaxWModel(pro)->getProblemTotalName() + "_numSamples_" + std::to_string(solIds.size()) + "_neighborK_" + std::to_string(neighborK);

	{
		std::ofstream nbnOut(dir + filename + "_nbn.txt");// (globalInfo.m_save_dir + filepath + "_nbn.txt");
		ofec::outputNBNdata(nbnOut, solIds, belong, belongDis, fitness);
		nbnOut.close();

	}


	{

		std::string filepath = dir;


		std::vector<TreeGraphSimple::NBNdata> nbn_data;
		transferNBNData(nbn_data, solIds, belong, belongDis, fitness);
		TreeGraphSimple nbn_graph;
		nbn_graph.setNBNdata(nbn_data);
		nbn_graph.calNetwork();



		std::string savepath = filepath + filename + "_network.txt";
		std::ofstream out(savepath);
		nbn_graph.outputNBNnetwork(out);
		out.close();


		ouputNBNdata(filepath + filename + "_nbn.txt", nbn_graph.get_nbn_data());

	}
}






void generatePro(std::shared_ptr<ofec::Environment>& env, ThreadInfoNBN& curcurInfo) {
	using namespace ofec;



	ParameterMap params;
	params["problem name"] = std::string("OneMaxWModel");
	params["number of variables"] = curcurInfo.m_dim;



	env.reset(Environment::create());
	env->recordInputParameters();
	env->initialize();

	env->setProblem(OneMaxWModel::create());
	env->problem()->inputParameters().input(params);
	//env->problem()->inputParameters().at("number of variables")->setValue(curcurInfo.m_dim);
	env->problem()->recordInputParameters();


	auto pro = env->problem();
	auto onemaxPro = GET_OneMaxWModel(pro);


	auto rnd = std::shared_ptr<Random>(new Random(0.5));;


	GET_OneMaxWModel(pro)->setInstanceId(curcurInfo.m_instance_id);
	GET_OneMaxWModel(pro)->setNeutralityMu(curcurInfo.m_neutrality_mu);
	GET_OneMaxWModel(pro)->setDummySelectRate(curcurInfo.m_dummy_select_rate);
	GET_OneMaxWModel(pro)->setEpistasisBlockSize(curcurInfo.m_epistasis_block_size);
	GET_OneMaxWModel(pro)->setRuggednessGamma(curcurInfo.m_ruggedness_gamma);



	pro->initialize(env.get());

}



void generateOneMaxProNeighborTask() {
	//std::vector<int> dims = { 80};
	//std::vector<int> ruggedness[2] = { {0,1, 2,4,6,8 ,10 } ,{0,1,4} };
	//std::vector<int> epistasis=  {0,1,2,3,4,5,6,7,8,9,10,20,31} ,{0,1,2,3,4} };
	//std::vector<int> neutral[2] = { { 0,1,2,3,5,7,9,10,20,30 },{0,1,4} };

	//	std::vector<std::vector<int>>  totalIds = { {1,0,0}, };

	using namespace ofec;

	ThreadGlobalNBN curGlobalInfo;
	ThreadInfoNBN curcurInfo;


	curGlobalInfo.m_save_dir = "//172.24.24.151/e/DiaoYiya/paper_com_experiment_data/oneMax/onemaxNeighborData_0531/";
	//	curGlobalInfo.m_save_dir = "/mnt/Data/Student/2018/YiyaDiao/NBN_data_paper_com/oneMax/onemaxNeighborData/";



	std::string dir = curGlobalInfo.m_save_dir;
	std::filesystem::create_directories(curGlobalInfo.m_save_dir);
	std::ofstream error_out(curGlobalInfo.m_save_dir + "_error_info.txt");
	curGlobalInfo.m_out.open(curGlobalInfo.m_save_dir + "finalResult.txt");
	curGlobalInfo.outputHead();
	curGlobalInfo.m_out << "ruggness" << "\t" << "neutrality" << "\t";
	curGlobalInfo.m_out << std::endl;


	int max_neutraility_mu = 60;
	int maxRugLoop = 5;
	int maxNeuLoop = 5;
	int maxELoop = 5;

	int maxDump = 20;

	std::vector<int> dims = { 120 };

	curcurInfo.m_instance_id = 1;
	curcurInfo.m_dim = 80;
	curcurInfo.m_ruggedness_gamma = 0;
	curcurInfo.m_epistasis_block_size = 0;
	curcurInfo.m_neutrality_mu = 0;
	curcurInfo.m_rug_ratio = 0;




	std::vector<ThreadInfoNBN> totalTasks;



	int taskId = 0;


	for (int curdim : dims) {

		curcurInfo.m_dim = curdim;
		curcurInfo.m_ruggedness_gamma = 0;
		curcurInfo.m_rug_ratio = 0;
		curcurInfo.m_neutrality_mu = 0;
		curcurInfo.m_epistasis_block_size = 0;

		for (int rugPar = 0; rugPar <= maxRugLoop; ++rugPar) {

			double ratio = double(rugPar) / double(maxRugLoop);
			curcurInfo.m_rug_ratio = ratio;
			curcurInfo.m_ruggedness_gamma = curcurInfo.m_dim * (curcurInfo.m_dim + 1) / 2.0 * ratio;
			totalTasks.push_back(curcurInfo);

		}


		curcurInfo.m_ruggedness_gamma = 0;
		curcurInfo.m_rug_ratio = 0;
		curcurInfo.m_neutrality_mu = 0;
		curcurInfo.m_epistasis_block_size = 0;

		for (int rugPar = 0; rugPar <= maxNeuLoop; ++rugPar) {

			double ratio = double(rugPar) / double(maxNeuLoop);

			curcurInfo.m_neutrality_mu = rugPar * max_neutraility_mu / maxNeuLoop;
			totalTasks.push_back(curcurInfo);
		}



		curcurInfo.m_ruggedness_gamma = 0;
		curcurInfo.m_rug_ratio = 0;
		curcurInfo.m_neutrality_mu = 0;
		curcurInfo.m_epistasis_block_size = 0;

		for (int rugPar = 0; rugPar <= maxELoop; ++rugPar) {


			curcurInfo.m_epistasis_block_size = rugPar * 4 + 2;
			totalTasks.push_back(curcurInfo);
		}

	}





	std::shared_ptr<ofec::Environment> env;
	std::shared_ptr<ofec::Random> rnd;

	int numSols(1e6);
	std::vector<std::shared_ptr<ofec::SolutionBase>> sols;

	for (auto& taskInfo : totalTasks) {
		generatePro(env, taskInfo);

		auto onemaxPro = CAST_ONEMAX(env->problem());

		std::cout << "running " << GET_OneMaxWModel(env->problem())->getProblemTotalName() << std::endl;
		rnd.reset(new ofec::Random(0.5));
		generateSolsRandom(env.get(), rnd.get(), sols, numSols);


		calNBNonexMaxNeighbor(env.get(), rnd.get(),
			dir, onemaxPro->numberVariables(), sols);

		int numNeighbor = onemaxPro->numberVariables();

		std::shared_ptr<ofec::SolutionBase> centerSol = sols.front();
		for (auto& it : sols) {
			if (centerSol->fitness() < it->fitness()) {
				centerSol = it;
			}
		}

		while (true) {
			numNeighbor /= 2.0;
			if (numNeighbor <= 3) {
				break;
			}

			std::cout << "neighbor k \t" << numNeighbor << std::endl;
			generateSols(env.get(), rnd.get(), *centerSol, sols, numNeighbor, numSols);


			calNBNonexMaxNeighbor(env.get(), rnd.get(),
				dir, numNeighbor, sols);
		}

	}

}



#include <filesystem>
#include <string>
#include <vector>
#include <iostream>

namespace ofec {

	void registerParamAbbr();
	void customizeFileName();
	void run() {
		registerParamAbbr();
		customizeFileName();
		using namespace std;


		generateOneMaxProNeighborTask();





		//runOneMaxTotalTaskAnalysis();

		//std::vector<std::string> filenames;
		//std::string dir = "//172.29.65.56/share/Student/2018/YiyaDiao/NBN_data/visualization_nbn_3/";

		//readDir(dir, filenames);

	//	calMultiTask();


		//calMultiTask();
		//analyseProTask();

		//filterMMOP_flacco();

		//filter_flacco_FreePeakBasinAtraction();

	}

	void registerParamAbbr() {}

	void customizeFileName() {}
}
