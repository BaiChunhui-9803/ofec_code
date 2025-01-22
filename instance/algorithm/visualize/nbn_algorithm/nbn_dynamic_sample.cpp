#include "nbn_dynamic_sample.h"
#include <chrono>
#include "../../../core/problem/continuous/continuous.h"
#include "../../../utility/nbn_visualization/nbn_edge_multithread_division.h"
//#include "sample2d_graph_algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#include <ui/custom/buffer/general/nearest_better_network/item_nearest_better_network.h>
//#include <ui/custom/buffer/algorithm/others/nbn_visualize/buffer_alg_nearest_better_network.h>
#endif
void ofec::NBN_DynamicSampleAlgorithm::sampleTwoObj()
{

	//std::vector<double> toParentNorDis(m_nbn_network.Sample().size(), 1e9);
	//std::vector<double> avgSampleNorDis(toParentNorDis.size(), 1e9);
	//std::pair<double, double> parDisRange(1e9, -1e9);
	//std::pair<double, double> avgDisRange(1e9, -1e9);

	//for (auto& it : m_nbn_network.Sample()) {

	//	if (!it->flagNoParent()) {
	//		parDisRange.first = std::min(it->m_disToParent, parDisRange.first);
	//		parDisRange.second = std::max(it->m_disToParent, parDisRange.second);

	//	}

	//	avgDisRange.first = std::min(it->m_neighbors_sum_dis, avgDisRange.first);
	//	avgDisRange.second = std::max(it->m_neighbors_sum_dis, avgDisRange.second);

	//}

	//for (int idx(0); idx < m_nbn_network.getSamples().size(); ++idx) {
	//	toParentNorDis[idx] = m_nbn_network.getSamples()[idx]->m_disToParent;
	//	avgSampleNorDis[idx] = m_nbn_network.getSamples()[idx]->m_neighbors_sum_dis;

	//	toParentNorDis[idx] = mapRealInside<double>(toParentNorDis[idx], parDisRange.first, parDisRange.second, 0, 1);
	//	avgSampleNorDis[idx] = mapRealInside<double>(avgSampleNorDis[idx], avgDisRange.first, avgDisRange.second, 0, 1);

	//}

	//double maxValue(0);
	//int sampleIdx(-1);
	//// function {a^3*b^3}
	//double curValue(0);
	//for (int idx(0); idx < avgSampleNorDis.size(); ++idx) {
	//	curValue = std::pow(toParentNorDis[idx], 3) * std::pow(avgSampleNorDis[idx], 3);
	//	if (maxValue < curValue) {
	//		maxValue = curValue;
	//		sampleIdx = idx;
	//	}
	//}


	//int num_vars = m_problem->numberVariables();
	//int number_objectives = m_problem->numberObjectives();
	//int num_cons = m_problem->numberConstraints();
	//m_cur_sol.reset(new SolutionType(number_objectives, num_cons, num_vars));

	//if (m_nbn_network.getSamples()[sampleIdx]->flagNoParent()) {
	//	m_problem->initializeSolution(*m_cur_sol, *m_nbn_network.getSamples()[sampleIdx]->m_cur_sol, m_problem->maximumVariableDistance(), m_random.get());
	//}
	//else {
	//	m_problem->initializeSolution(*m_cur_sol, *m_nbn_network.getSamples()[sampleIdx]->m_cur_sol, m_nbn_network.getSamples()[sampleIdx]->m_disToParent, m_random.get());
	//}

	//m_eval_fun(*m_cur_sol, m_problem.get());
	//m_nbn_network.addRandomSol(m_cur_sol);

}

void ofec::NBN_DynamicSampleAlgorithm::outputNBNdata() {

	
	std::vector<int> markerPonits;
	std::vector<std::shared_ptr<ofec::SolutionBase>> samples;
	std::vector<int> belongs;
	std::vector<double> fitness;
	std::vector<double> dis2parent;
	std::vector<bool> optFlag;

	auto& nbn_division(m_nbn_network.m_division);
	nbn_division->getSharedNBN(samples, fitness, belongs, dis2parent, optFlag);

	//if (samples.empty())return;

	auto& nbn_data(m_nbn_data);
	nbn_data.initializeBasicInfo(m_problem.get());
	nbn_data.update(samples, belongs, fitness, optFlag);
	nbn_data.updateFileInfo(*m_param, *m_param);

	nbn_data.m_nbn_data.m_abfilepath = m_fileDir + "/" + nbn_data.m_nbn_data.m_filename + ".txt";

	nbn_data.m_nbn_data.output(nbn_data.m_nbn_data.m_abfilepath);
}


void ofec::NBN_DynamicSampleAlgorithm::outputNBNdata(const std::string& path_dir) {
	std::vector<int> markerPonits;
	std::vector<std::shared_ptr<ofec::SolutionBase>> samples;
	std::vector<int> belongs;
	std::vector<double> fitness;
	std::vector<double> dis2parent;
	std::vector<bool> optFlag;

	auto& nbn_division(m_nbn_network.m_division);
	nbn_division->getSharedNBN(samples, fitness, belongs, dis2parent, optFlag);

	//if (samples.empty())return;

	auto& nbn_data(m_nbn_data);
	nbn_data.initializeBasicInfo(m_problem.get());
	nbn_data.update(samples, belongs, fitness, optFlag);
	nbn_data.updateFileInfo(*m_param, *m_param);
	nbn_data.m_nbn_data.m_abfilepath = path_dir + "/" + nbn_data.m_nbn_data.m_filename + ".txt";
	
	nbn_data.m_nbn_data.output(nbn_data.m_nbn_data.m_abfilepath);
}
	
void ofec::NBN_DynamicSampleAlgorithm::sample1() {
	int num_vars = m_problem->numberVariables();
	int number_objectives = m_problem->numberObjectives();
	int num_cons = m_problem->numberConstraints();

	m_cur_sol.reset(m_problem->createSolution());
	//m_cur_sol.reset(new SolutionType(number_objectives, num_cons, num_vars));
	m_problem->initializeSolution(*m_cur_sol, m_random.get());
	m_eval_fun(*m_cur_sol, m_problem.get());
	int id(0);
	m_nbn_network.m_division->addSol(*m_cur_sol, id);

//	m_nbn_network.addRandomSol(m_cur_sol,id);
}

void ofec::NBN_DynamicSampleAlgorithm::samppleCEC2015FO5() {
	int num_vars = m_problem->numberVariables();
	int number_objectives = m_problem->numberObjectives();
	int num_cons = m_problem->numberConstraints();

	m_cur_sol.reset(m_problem->createSolution());
	//m_cur_sol.reset(new SolutionType(number_objectives, num_cons, num_vars));
	m_problem->initializeSolution(*m_cur_sol, m_random.get());

	auto cur_sol(dynamic_cast<Solution<>*>(m_cur_sol.get()));

	// sample for cec2015f05 

	if (m_problem->name() == "MMOP_CEC2015_F05" && m_problem->numberVariables() == 2)
	{

		//	cout << "sample of MMOP_CEC2015_F05" << std::endl;
		if (m_random->uniform.next() < 0.5) {
			double ratio = m_random->uniform.next();
			cur_sol->variable()[0] = ratio * 200 - 100;
			cur_sol->variable()[1] = -30 - ratio * 30 - m_random->uniform.nextNonStd<double>(-30, 0);
		}
		else {
			double ratio = m_random->uniform.next();
			cur_sol->variable()[0] = -45 + 30 * ratio + m_random->uniform.nextNonStd<double>(0, 30);
			cur_sol->variable()[1] = ratio * 200 - 100;

		}
	}


	m_eval_fun(*m_cur_sol, m_problem.get());
	int sol_id(0);
//	m_nbn_network.addRandomSol(m_cur_sol,sol_id);


	
	m_nbn_network.m_division->addSol(*m_cur_sol, sol_id);

}


void ofec::NBN_DynamicSampleAlgorithm::initialize_()
{


	Algorithm::initialize_();
	m_eval_fun =
		[](SolutionBase& sol, Problem *pro) {
		using namespace ofec;
		sol.evaluate(pro, -1, false);
		ofec::Real pos = (pro->optimizeMode(0) == ofec::OptimizeMode::kMaximize) ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
	};

	auto& v = *m_param;

	m_continousSampleSize = v.get<int>("continous sample size", 1e5);
	m_initSampleSize = v.get<int>("init sample size", mc_initSampleSize);
	m_filter_flag = v.get<bool>("NBN filter", false);
	m_opt_flag = v.get<bool>("add optimal", false);
	
	//m_continousSampleSize = v.has("continous sample size") ? v.get<int>("continous sample size") : 1e5;
	//m_initSampleSize = v.has("init sample size") ? v.get<int>("init sample size") : mc_initSampleSize;

	//m_filter_flag = v.has("NBN filter") ? v.get<bool>("NBN filter") : false;
	//m_opt_flag = v.has("add optimal") ? v.get<bool>("add optimal") : false;



	//UTILITY::generateSamplesEval<Solution<>, std::unique_ptr>(
	//	sols, rnd, pro, eval_fun);
	int init_sample_num(0);
	bool m_flag_grid = false;
	if (m_problem->hasTag(ofec::ProblemTag::kConOP)) {
		init_sample_num = m_continousSampleSize;
		m_flag_grid = true;
	}
	else {
		init_sample_num = m_initSampleSize;
		//initSolsRand();
	}
//	m_nbn_network.initialize(m_random.get(), m_problem.get(),init_sample_num,m_eval_fun, m_flag_grid);
	
	//if (m_opt_flag) {
	//	addOpts();
	//}

#ifdef OFEC_DEMO
	ofec::NBN_DivisionData::initialize_static_info();
	ofec::NBN_Data::updateDirPath();
#endif // DEBUG


	m_startTime = std::chrono::system_clock::now();
//	Random *rnd = ADD_RND(0.3);
	m_nbn_network.initialize_division(m_problem.get(), m_random.get(), v);

	//std::cout << "init sample size\t" << m_nbn_network.m_division->size() << std::endl;


#ifdef  OFEC_DEMO
	ofec::NBN_DivisionData::ms_division = m_nbn_network.m_division;
#endif
//	

	updateBuffer();
	//updateBuffer();
	outputNBNdata();

}

void ofec::NBN_DynamicSampleAlgorithm::addOpts() {

//	m_nbn_network.addOpt();

//	auto& optBase(m_problem->optBase());
//	std::cout << "dynamic sample alg opt numbers\t" << optBase.numberVariables() << std::endl;

}

//
//void ofec::NBN_DynamicSampleAlgorithm::initSolsRand() {
//	//nbn.setEvalFun(eval_fun);
//	m_nbn_network.initRandomSol(m_initSampleSize, m_eval_fun);
//	updateBuffer();
//}
//
//
//
//void ofec::NBN_DynamicSampleAlgorithm::initSols2Dgraph() {
//
//
//	Sample2D_Graph_Algorithm sample_alg;
//	sample_alg.setMaxSampleSize(m_continousSampleSize);
//	sample_alg.init(m_problem.get(), m_random.get(), m_eval_fun);
//	sample_alg.calculate();
//	if (m_opt_flag) {
//		sample_alg.insertOpt();
//	}
//
//	sample_alg.udpate_network();
//
//#ifdef  OFEC_DEMO
//	//m_nbn_network.outputData(ofec_demo::bufferAlgNearestBetterNetwork::ms_NBN_record);
//
//	auto& samples(ofec_demo::bufferAlgNearestBetterNetwork::ms_samples);
//	auto& belongs(ofec_demo::bufferAlgNearestBetterNetwork::ms_belongs);
//	auto& fitness(ofec_demo::bufferAlgNearestBetterNetwork::ms_fitness);
//
//	auto& neighbors(ofec_demo::bufferAlgNearestBetterNetwork::ms_neighbors);
//
//
//	if (m_filter_flag) {
//		sample_alg.getNearestBetterNetworkShareMemoryFilter(samples, fitness, belongs);
//	}
//	else {
//		sample_alg.getNearestBetterNetworkShareMemory(samples, fitness, belongs);
//	}
//
//	neighbors = sample_alg.neighbors();
//
//	ofec_demo::g_buffer->appendAlgBuffer(this);
//	//	updateBuffer();
//
//#endif //  OFEC_DEMO
//
//
//
//	//SampleNearestBetterNetworkRecord record;
//	//record.m_problem.get() = m_nbn_network.idPro();
//	//record.m_random.get() = m_nbn_network.idRnd();
//	//record.m_com_fun = m_nbn_network.comFun();
//	//sample_alg.outputData(record);
//	//m_nbn_network.insertData(record);
//
//
//
//}

void ofec::NBN_DynamicSampleAlgorithm::run_()
{

	auto from_time = std::chrono::system_clock::now();
	auto to_time = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(to_time - from_time);
	double duration_seconds = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

	std::shared_ptr<SolutionBase> cur_sol;
	//int num_vars = m_problem->numberVariables();
	//int number_objectives = m_problem->numberObjectives();
	//int num_cons = m_problem->numberConstraints();

	int numLoop(1);
	

	while (!isFinish()) {
		//if(dynamic_cast<ofec::NBN_EdgeMultiThreadDivision&>(*m_nbn_network.m_division))
		
		++numLoop;

		if (m_nbn_network.m_division_type == NBN_VisualizationData::NBN_Division_Type::kEdgeMultiThreadDivision) {
			dynamic_cast<ofec::NBN_EdgeMultiThreadDivision&>(*m_nbn_network.m_division).updateDivisionMultithread();

			std::cout << "number loop\t" << numLoop << std::endl;
			updateBuffer();
		}
		else if (m_nbn_network.m_division_type == NBN_VisualizationData::NBN_Division_Type::kEdgeDivision) {
			dynamic_cast<ofec::NBN_EdgeDivision&>(*m_nbn_network.m_division).updateDivision();

			std::cout << "number loop\t" << numLoop << std::endl;
			updateBuffer();
		}
		else {
			cur_sol.reset(m_problem->createSolution());
			m_problem->initializeSolution(*cur_sol, m_random.get());
			m_eval_fun(*cur_sol, m_problem.get());
			int belong_id(-1);
			m_nbn_network.m_division->addSol(*cur_sol,belong_id);

			to_time = std::chrono::system_clock::now();
			duration = std::chrono::duration_cast<std::chrono::microseconds>(to_time - from_time);
			duration_seconds = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

			if (duration_seconds > m_duration_seconds) {

				//	std::cout << "dynamic sample algorithm current points\t" << m_nbn_network.m_division->size() << std::endl;

				from_time = std::chrono::system_clock::now();

				std::cout << "number loop\t" << numLoop << std::endl;
				updateBuffer();
			}
		}
		

	
		//sample1();



		//samppleCEC2015FO5();



		//// sample for cec2015f05 

		//if(m_problem->name()=="MMOP_CEC2015_F05"&&m_problem->numberVariables()==2)
		//{

		////	cout << "sample of MMOP_CEC2015_F05" << std::endl;
		//	if (m_random->uniform.next() < 0.5) {
		//		double ratio = m_random->uniform.next();
		//		cur_sol->variable()[0] = ratio * 200 - 100;
		//		cur_sol->variable()[1] = -30 - ratio * 30 - m_random->uniform.nextNonStd<double>(-30,0);
		//	}
		//	else {
		//		double ratio = m_random->uniform.next();
		//		cur_sol->variable()[0] = -45+ 30*ratio+ m_random->uniform.nextNonStd<double>(0, 30);
		//		cur_sol->variable()[1] = ratio * 200 - 100;

		//	}
		//}


	



	}
}


void ofec::NBN_DynamicSampleAlgorithm::updateBuffer() {

#ifdef  OFEC_DEMO
	/*
	namespace ofec_demo {
	/* Buffer of Sample Algorithm with NBN */

	//auto& samples(ofec_demo::ItemNearestBetterNetwork::ms_samples);
	//auto& belongs(ofec_demo::ItemNearestBetterNetwork::ms_belongs);
	//auto& fitness(ofec_demo::ItemNearestBetterNetwork::ms_fitness);
	//auto& marker(ofec_demo::ItemNearestBetterNetwork::ms_markerPonits);
	//m_nbn_network.getSharedData(samples, fitness, belongs);


	//ofec_demo::bufferAlgNearestBetterNetwork::ms_NBN_record.getNearestBetterNetwork(samples, fitness, belongs);

	ofec_demo::g_buffer->appendAlgBuffer(this);

	

#endif //  OFEC_DEMO
}