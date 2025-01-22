#include"sample_nearest_better_network.h"
#include "../../../utility/function/custom_function.h"
#include "../../../core/problem/continuous/continuous.h"
#include <thread>

double ofec::NBNStaticPriorityQueue::avgDis() {
	if (empty())return 1e9;
	else {
		double avgDis(0);
		for (int idx(0); idx < m_cur_size; ++idx) {
			avgDis += m_cur_sol->m_cur_sol->variableDistance(*m_data[idx]->m_cur_sol, m_problem.get());
		}
		return avgDis / m_cur_size;
	}
}

bool ofec::SampleNearestBetterNetwork::addRandomSol(std::shared_ptr<SolutionType>& sol, bool flag_opt )
{
	for (auto& it : m_samples) {
		if (m_problem->same(*it->m_cur_sol, *sol)) {
			return false;
		}
	}
	auto& curNode(insertNode(sol));
	curNode->m_flag_opt = flag_opt;
	auto listIter(m_sorted_list.front());
	auto beforeIter(m_sorted_list.head());
	while (!m_sorted_list.isEnd(listIter)) {
		listIter->m_info->addNeighbor(m_problem.get(), curNode);
		curNode->addNeighbor(m_problem.get(), listIter->m_info);

		if (curNode->m_cur_sol->fitness() > 
			listIter->m_info->m_cur_sol->fitness()) {
			listIter->m_info->updateParent(curNode, m_problem.get());
		}
		//else if(curNODEIN)
		else if (curNode->m_cur_sol->fitness() ==
			listIter->m_info->m_cur_sol->fitness()) {
		}
		else {
			curNode->updateParent(listIter->m_info, m_problem.get());
			beforeIter = listIter;
		}
		listIter = listIter->m_after;
	}
	NearestBetterNetworkList::insertNode(beforeIter, &curNode->m_sorted_list_node);
	return true;
}

void ofec::SampleNearestBetterNetwork::addNearSol(const std::shared_ptr<NearestBetterNetworkNode>& center, std::unique_ptr<SolutionType>& sol)
{
}

void ofec::SampleNearestBetterNetwork::initRandomSol(int numSols, const std::function<void(SolutionType& sol, Problem *pro)>& eval_fun)
{
	std::vector<std::shared_ptr<SolutionType>> sols(numSols);
	sampleRandom(sols,eval_fun);
	for (int idx(0); idx < sols.size(); ++idx) {
		addRandomSol(sols[idx]);
	}
}

void ofec::SampleNearestBetterNetwork::sampleRandom(std::vector<std::shared_ptr<SolutionType>>& sols, const std::function<void(SolutionType& sol, Problem *pro)>& eval_fun)
{
	
	//UTILITY::generateSamplesEval<SolutionType, std::shared_ptr>
	//	(sols, GET_RND(m_random.get()), m_problem.get(), eval_fun);



	//int num_task = std::thread::hardware_concurrency();
	//int num_samples = sols.size();
	//std::vector<std::thread> thrds;
	//std::vector<int> tasks;
	//UTILITY::assignThreads(num_samples, num_task, tasks);
	//std::pair<int, int> from_to;
	//for (size_t i = 0; i < num_task; ++i) {
	//	from_to.first = tasks[i];
	//	from_to.second = tasks[i + 1];
	//	double randSeed(double(i + 1) / double(num_task + 1));
	//	thrds.push_back(std::thread(
	//		&SampleNearestBetterNetwork::sampleRandomThreadTask, this,
	//		std::ref(sols), from_to, randSeed));
	//}
	//for (auto& thrd : thrds)
	//	thrd.join();

}

void ofec::SampleNearestBetterNetwork::sampleRandomThreadTask(
	std::vector<std::shared_ptr<SolutionType>>& sols, 
	std::pair<int, int> from_idx, double randSeed)
{
	Random rand(randSeed);

	for (int idx(from_idx.first); idx < from_idx.second; ++idx) {
		sols[idx].reset(new SolutionType(m_problem->numberObjectives(),m_problem->numberConstraints(),m_problem->numberVariables()));
		m_problem->initializeSolution(*sols[idx], rnd);
	}
}

void ofec::SampleNearestBetterNetwork::getNearestBetterNetwork(
	std::vector<std::shared_ptr<SolutionType>>& sols,
	std::vector<double>& fitness,
	std::vector<int>& belong) {
	sols.resize(m_samples.size());
	belong.resize(m_samples.size());
	fitness.resize(m_samples.size());
	int cur_id(0);
	auto listIter(m_sorted_list.front());
	while (!m_sorted_list.isEnd(listIter)) {

		listIter->m_info->m_sorted_id = cur_id;
		sols[cur_id].reset(new SolutionType(*listIter->m_info->m_cur_sol.get()));
		if (listIter->m_info->m_parent != nullptr) {
			belong[cur_id] = cur_id;
		}
		else {
			belong[cur_id] = listIter->m_info->m_parent->m_sorted_id;
		}
		fitness[cur_id] = sols[cur_id]->fitness();
		++cur_id;
		listIter = listIter->m_after;
	}
}

void ofec::SampleNearestBetterNetwork::getNearestBetterNetwork(
	std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
	std::vector<double>& fitness,
	std::vector<int>& belong)
{
	sols.resize(m_samples.size());
	belong.resize(m_samples.size());
	fitness.resize(m_samples.size());
	int cur_id(0);
	auto listIter(m_sorted_list.front());
	while (!m_sorted_list.isEnd(listIter)) {
		
		listIter->m_info->m_sorted_id = cur_id;
		sols[cur_id].reset(new SolutionType(*listIter->m_info->m_cur_sol.get()));
		if (listIter->m_info->m_parent==nullptr) {
			belong[cur_id] = cur_id;
		}
		else {
			belong[cur_id] = listIter->m_info->m_parent->m_sorted_id;
		}
		fitness[cur_id] = sols[cur_id]->fitness();
		++cur_id;
		listIter = listIter->m_after;
	}
}
void ofec::SampleNearestBetterNetwork::getNearestBetterNetworkShareMemory(
	std::vector<std::shared_ptr<SolutionType>>& sols,
	std::vector<double>& fitness,
	std::vector<int>& belong) {
	sols.resize(m_samples.size());
	belong.resize(m_samples.size());
	fitness.resize(m_samples.size());
	int cur_id(0);
	auto listIter(m_sorted_list.front());
	while (!m_sorted_list.isEnd(listIter)) {

		listIter->m_info->m_sorted_id = cur_id;
		sols[cur_id] = listIter->m_info->m_cur_sol;
		if (listIter->m_info->m_parent == nullptr) {
			belong[cur_id] = cur_id;
		}
		else {
			belong[cur_id] = listIter->m_info->m_parent->m_sorted_id;
		}
		fitness[cur_id] = sols[cur_id]->fitness();
		++cur_id;
		listIter = listIter->m_after;
	}
}



void ofec::SampleNearestBetterNetwork::outputData(SampleNearestBetterNetworkRecord& record) {
	 record.m_random.get() = m_random.get();
	 record.m_problem.get() = m_problem.get();
	 record.m_com_fun = m_com_fun;
	// record.m_cur_visited_id = 0;

	record.m_samples.resize(m_samples.size());
	for (int idx(0); idx < m_samples.size(); ++idx) {
		record.m_samples[idx].reset(new SolutionType(*m_samples[idx]->m_cur_sol));
	}
	record.m_sorted_ids.clear();
	auto listIter(m_sorted_list.front());
	while (!m_sorted_list.isEnd(listIter)) {
		record.m_sorted_ids.push_back(listIter->m_info->m_sample_id);
		listIter = listIter->m_after;
	}

	record.m_parent_ids.resize(m_samples.size());
	record.m_neighbor_ids.resize(m_samples.size());
	for (int idx(0); idx < m_samples.size(); ++idx) {
		auto& curNode(m_samples[idx]);
		auto& curNeiIdx(record.m_neighbor_ids[idx]);
		for (int neiIdx(0); neiIdx < curNode->m_neighbors.size(); ++neiIdx) {
			auto& neiNode(curNode->m_neighbors.data()[neiIdx]);
			curNeiIdx.push_back(neiNode->m_sample_id);
		}
		if (curNode->m_parent == nullptr) {
			record.m_parent_ids[idx] = idx;
		}
		else {
			record.m_parent_ids[idx] = curNode->m_parent->m_sample_id;
		}
	//	record.m_parent_ids[idx] = curNode->m_parent->m_sample_id;
	}
	//record

}
void ofec::SampleNearestBetterNetwork::insertData(const SampleNearestBetterNetworkRecord& record) {
	clear();
	m_random.get() = record.m_random.get();
	m_problem.get() = record.m_problem.get();
	m_com_fun = record.m_com_fun;
	m_cur_visited_id = 0;

	std::shared_ptr<SolutionType> cur_sol;

	for (int idx(0); idx < record.m_samples.size(); ++idx)
	{
		cur_sol.reset(new SolutionType(*record.m_samples[idx]));
		insertNode(cur_sol);
	}

	auto listIter(m_sorted_list.head());
	for (int idx(0); idx < record.m_samples.size(); ++idx) {
		int cur_id = record.m_sorted_ids[idx];
		auto& curNode(m_samples[cur_id]);
		NearestBetterNetworkList::insertNode(listIter, &curNode->m_sorted_list_node);
		listIter = listIter->m_after;
		//std::cout << "cur_id\t" << cur_id << std::endl;
	}
	for (int idx(0); idx < record.m_samples.size(); ++idx) {
		auto& curNode(m_samples[idx]);
		auto& parentNode(m_samples[record.m_parent_ids[idx]]);
		if(curNode->m_sample_id!=parentNode->m_sample_id)
		curNode->updateParent(parentNode, m_problem.get());
		for (auto& neighbor_idx : record.m_neighbor_ids[idx]) {
			auto& neiNode(m_samples[neighbor_idx]);
			curNode->addNeighbor(m_problem.get(), neiNode);
		}
	}
}
