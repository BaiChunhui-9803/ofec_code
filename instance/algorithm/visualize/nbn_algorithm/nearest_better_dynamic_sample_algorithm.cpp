#include "nearest_better_dynamic_sample_algorithm.h"

#include <chrono>
#include "../../../core/problem/continuous/continuous.h"
#include "sample2d_graph_algorithm.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#include <ui/custom/buffer/algorithm/others/nbn_visualize/buffer_alg_nearest_better_network.h>
#endif
//void ofec::NearestBetterDynamicSampleAlgorithm::record()
//{
//
//}

void ofec::NearestBetterDynamicSampleAlgorithm::sampleTwoObj()
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


void ofec::NearestBetterDynamicSampleAlgorithm::sample1() {
	int num_vars = m_problem->numberVariables();
	int number_objectives = m_problem->numberObjectives();
	int num_cons = m_problem->numberConstraints();

	m_cur_sol.reset(m_problem->createSolution());
	//m_cur_sol.reset(new SolutionType(number_objectives, num_cons, num_vars));
	m_problem->initializeSolution(*m_cur_sol, m_random.get());

	m_eval_fun(*m_cur_sol, m_problem.get());
	m_nbn_network.addRandomSol(m_cur_sol);
}

void ofec::NearestBetterDynamicSampleAlgorithm::samppleCEC2015FO5() {
	int num_vars = m_problem->numberVariables();
	int number_objectives = m_problem->numberObjectives();
	int num_cons = m_problem->numberConstraints();

	m_cur_sol.reset(m_problem->createSolution());
	//m_cur_sol.reset(new SolutionType(number_objectives, num_cons, num_vars));
	m_problem->initializeSolution(*m_cur_sol, m_random.get());

	auto cur_sol(dynamic_cast<Solution<>* >(m_cur_sol.get()));

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
	m_nbn_network.addRandomSol(m_cur_sol);
}


void ofec::NearestBetterDynamicSampleAlgorithm::initialize_()
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
	m_continousSampleSize = v.has("continous sample size") ? v.get<int>("continous sample size") : 1e5;
	m_initSampleSize = v.has("init sample size") ? v.get<int>("init sample size") : mc_initSampleSize;

	m_filter_flag = v.has("NBN filter") ? v.get<bool>("NBN filter") : false;
	m_opt_flag = v.has("add optimal") ? v.get<bool>("add optimal") : false;


	if (v.has("NBN filter")) {
		std::cout << "NBN filter" << std::endl;
	}

	if (v.has("add optimal")) {
		std::cout << "add optimal" << std::endl;
	}


	//UTILITY::generateSamplesEval<Solution<>, std::unique_ptr>(
	//	sols, rnd, pro, eval_fun);
	m_nbn_network.initialize(m_random.get(), m_problem.get());
	if (m_opt_flag) {
		addOpts();
	}

	initSolsRand();

	//if (m_problem->hasTag(ofec::ProblemTag::kConOP)) {
	//	initSols2Dgraph();
	//}
	//else {
	//	initSolsRand();
	//}

	//updateBuffer();


}

void ofec::NearestBetterDynamicSampleAlgorithm::addOpts() {

	m_nbn_network.addOptSols(m_eval_fun);

	//std::shared_ptr<SolutionType> cur_sol;
	//int num_vars = m_problem->numberVariables();
	//int number_objectives = m_problem->numberObjectives();
	//int num_cons = m_problem->numberConstraints();

	auto& optBase(m_problem->optBase());
	//std::vector<std::shared_ptr<SolutionType>> opt_sols;
	//for (size_t idx(0); idx < optBase.numberVariables(); ++idx) {
	//	cur_sol.reset(new SolutionType(number_objectives, num_cons, num_vars));
	//	//cur_sol->variable() = dynamic_cast<const SolutionType&>(optBase.variable(idx)).variable();
	//	m_eval_fun(*cur_sol, m_problem.get());
	//	opt_sols.push_back(cur_sol);
	//	m_nbn_network.addRandomSol(cur_sol);

	//}
	std::cout << "dynamic sample alg opt numbers\t" << optBase.numberVariables() << std::endl;

//#ifdef  OFEC_DEMO
//
//	auto& markderIdxs(ofec_demo::bufferAlgNearestBetterNetwork::ms_markerPonits);
//	markderIdxs.clear();
//	for (size_t idx(0); idx < optBase.numberVariables(); ++idx) {
//		markderIdxs.push_back(idx);
//	}
//
//#endif //  OFEC_DEMO
}


void ofec::NearestBetterDynamicSampleAlgorithm::initSolsRand() {
	//nbn.setEvalFun(eval_fun);
	m_nbn_network.initRandomSol(m_initSampleSize, m_eval_fun); 
	updateBuffer();
}



void ofec::NearestBetterDynamicSampleAlgorithm::initSols2Dgraph() {


	Sample2D_Graph_Algorithm sample_alg;
	sample_alg.setMaxSampleSize(m_continousSampleSize);
	sample_alg.init(m_problem.get(), m_random.get(), m_eval_fun);
	sample_alg.calculate();
	if (m_opt_flag) {
		sample_alg.insertOpt();
	}

	sample_alg.udpate_network();

#ifdef  OFEC_DEMO
	//m_nbn_network.outputData(ofec_demo::bufferAlgNearestBetterNetwork::ms_NBN_record);

	auto& samples(ofec_demo::bufferAlgNearestBetterNetwork::ms_samples);
	auto& belongs(ofec_demo::bufferAlgNearestBetterNetwork::ms_belongs);
	auto& fitness(ofec_demo::bufferAlgNearestBetterNetwork::ms_fitness);

	auto& neighbors(ofec_demo::bufferAlgNearestBetterNetwork::ms_neighbors);


	if (m_filter_flag) {
		sample_alg.getNearestBetterNetworkShareMemoryFilter(samples, fitness, belongs);
	}
	else {
		sample_alg.getNearestBetterNetworkShareMemory(samples, fitness, belongs);
	}

	neighbors = sample_alg.neighbors();

	ofec_demo::g_buffer->appendAlgBuffer(this);
//	updateBuffer();

#endif //  OFEC_DEMO



	//SampleNearestBetterNetworkRecord record;
	//record.m_problem.get() = m_nbn_network.idPro();
	//record.m_random.get() = m_nbn_network.idRnd();
	//record.m_com_fun = m_nbn_network.comFun();
	//sample_alg.outputData(record);
	//m_nbn_network.insertData(record);



}

void ofec::NearestBetterDynamicSampleAlgorithm::run_()
{

	auto from_time = std::chrono::system_clock::now();
	auto to_time = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(to_time - from_time);
	double duration_seconds = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

	std::shared_ptr<SolutionBase> cur_sol;
	//int num_vars = m_problem->numberVariables();
	//int number_objectives = m_problem->numberObjectives();
	//int num_cons = m_problem->numberConstraints();


	while (!terminating()) {

		sample1();
		//samppleCEC2015FO5();

		//cur_sol.reset(new SolutionType(number_objectives, num_cons, num_vars));
		//m_problem->initializeSolution(*cur_sol, m_random.get());

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


		//m_eval_fun(*cur_sol, m_problem.get());
		//m_nbn_network.addRandomSol(cur_sol);

		to_time = std::chrono::system_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(to_time - from_time);
		duration_seconds = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

		if (duration_seconds > m_duration_seconds) {

			std::cout << "dynamic sample algorithm current points\t" << m_nbn_network.Sample().size() << std::endl;

			from_time = std::chrono::system_clock::now();

			updateBuffer();
		}

	}
}


void ofec::NearestBetterDynamicSampleAlgorithm::updateBuffer() {

#ifdef  OFEC_DEMO
	auto& samples(ofec_demo::bufferAlgNearestBetterNetwork::ms_samples);
	auto& belongs(ofec_demo::bufferAlgNearestBetterNetwork::ms_belongs);
	auto& fitness(ofec_demo::bufferAlgNearestBetterNetwork::ms_fitness);
	auto& marker(ofec_demo::bufferAlgNearestBetterNetwork::ms_markerPonits);
	m_nbn_network.getNearestBetterNetwork(samples, fitness, belongs, marker);


	//ofec_demo::bufferAlgNearestBetterNetwork::ms_NBN_record.getNearestBetterNetwork(samples, fitness, belongs);

	ofec_demo::g_buffer->appendAlgBuffer(this);

#endif //  OFEC_DEMO
}