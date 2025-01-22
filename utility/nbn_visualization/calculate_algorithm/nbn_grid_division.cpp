#include "nbn_grid_division.h"
#include "../../core/problem/continuous/continuous.h"
#include "../function/custom_function.h"



void ofec::NBN_GridDivision::initialize_(bool flag_grid_sample )  {
	m_dim = m_pro->numVariables();
	initNodesNormal();
	updateFitness();
	calculate();
}


void ofec::NBN_GridDivision::idxToVec(int idx, std::vector<int>& cur) {

	for (int idDim(m_dim - 1); idDim >= 0; --idDim) {
		cur[idDim] = idx % m_dim_div[idDim];
		idx /= m_dim_div[idDim];
	}

}
void ofec::NBN_GridDivision::vecToIdx(const std::vector<int>& cur, int& idx) {
	//for(int idDim)
	idx = 0;
	for (int idDim(0); idDim < m_dim; ++idDim) {
		idx *= m_dim_div[idDim];
		idx += cur[idDim];
	}

}


void ofec::NBN_GridDivision::solToVec(const SolBase& sol, std::vector<int>& vec) {
	auto& cur_sol(dynamic_cast<const Solution<>&>(sol));
	vec.resize(m_dim);

	for (int idDim(0); idDim < m_dim; ++idDim) {
		vec[idDim] = (cur_sol.variable()[idDim] - m_range_from[idDim]) / m_range_div[idDim];
	//	cur_sol->variable()[idDim] = boundary[idDim].first + curVec[idDim] * range_div[idDim];
	}
}
bool ofec::NBN_GridDivision::judgeFeasible(const std::vector<int>& cur) {
	for (int idx(0); idx < m_dim; ++idx) {
		if (cur[idx] >= 0 && cur[idx] < m_dim_div[idx]) {

		}
		else return false;
	}
	return true;
} 

void ofec::NBN_GridDivision::initNodesNormal() {

	//std::vector<Node*> nodes;

	int numSample_int = m_maxSample;
	m_dim_div.resize(m_dim);
	int dim_div = std::max<int>(2,exp(log(double(numSample_int)) / double(m_dim)));
	numSample_int = 1;
	for (auto& it : m_dim_div) {
		it = dim_div;
		numSample_int *= it;
	}

	m_nbn_nodes.resize(numSample_int);
	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
		m_nbn_nodes[idx].init(idx);
		auto& cur_sol(m_nbn_nodes[idx].m_sol);
		cur_sol.reset(new Solution<>(
			m_pro->numObjectives(),
			m_pro->numConstraints(),
			m_pro->numVariables()));
//		cur_sol->setId(idx);
	//	m_nbn_nodes[idx].m_id = idx;
	}


	//for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
	//	m_nbn_nodes[idx].m_sol = m_cur_sols[m_nbn_nodes[idx].m_sol_id];
	//}

	std::vector<double> range_div(m_dim);
	auto boundary = CAST_CONOP(m_pro.get())->boundary();
	for (int idDim(0); idDim < m_dim; ++idDim) {
		auto& range_var = CAST_CONOP(m_pro.get())->range(idDim);
		range_div[idDim] = (range_var.second - range_var.first) / double(m_dim_div[idDim]);
	}

	m_range_div = range_div;
	m_range_from.resize(range_div.size());
	for (int idDim(0); idDim < range_div.size(); ++idDim) {
		auto& range_var = CAST_CONOP(m_pro.get())->range(idDim);
		m_range_from[idDim] = range_var.first;
	}


	std::vector<int> curVec(m_dim, 0);
	int curIdx(0);
	// set solutions
	for (size_t idx(0); idx < numSample_int; ++idx) {
		idxToVec(idx, curVec);
		vecToIdx(curVec, curIdx);
		if (idx != curIdx) {
			THROW("error at wrong vec to idx \n");
			//std::cout << "error" << endl;
		}
		m_nbn_nodes[idx].m_vec = curVec;
		auto& cur_sol(dynamic_cast<Solution<>&>(*m_nbn_nodes[idx].m_sol));
		for (int idDim(0); idDim < m_dim; ++idDim) {
			cur_sol.variable()[idDim] = boundary[idDim].first + curVec[idDim] * range_div[idDim];
		}
	}

	// set neighbor sols
	std::vector<int> neiVec(m_dim, 0);
	int neiIdx(0);
	for (int idx(0); idx < numSample_int; ++idx) {
		auto& cur_node = m_nbn_nodes[idx];
		idxToVec(idx, curVec);
		curIdx = idx;
		neiVec = curVec;
		for (int idDim(0); idDim < m_dim; ++idDim) {
			if (--neiVec[idDim] >= 0) {
				vecToIdx(neiVec, neiIdx);
				auto& nei_node = m_nbn_nodes[neiIdx];
				double cur_dis = m_pro->norVariableDistance(*cur_node.m_sol, *nei_node.m_sol);
				cur_node.m_neighbor_id_dis.emplace_back(neiIdx, cur_dis);
				nei_node.m_neighbor_id_dis.emplace_back(curIdx, cur_dis);
				//std::pair<int, double> nei_info = {neiIdx.}
				//cur_node.m_neighbor.push_back(neiIdx);
				//cur_node.m_neighbor_dis.push_back(cur_dis);
				//nei_node.m_neighbor.push_back(curIdx);
				//nei_node.m_neighbor_dis.push_back(cur_dis);
				m_neighbors.push_back({ idx, neiIdx });
			}
			++neiVec[idDim];
		}
	}

	if (m_pro->hasTag(ofec::ProTag::kConOP) && m_pro->numVariables() == 2) {
		m_fitness_landscape_idxs.front().resize(m_dim_div.front());
		m_fitness_landscape_idxs.back().resize(m_dim_div.back());
		for (size_t idx(0); idx < numSample_int; ++idx) {
			idxToVec(idx, curVec);
			m_fitness_landscape_idxs[0][curVec[0]] = idx;
			m_fitness_landscape_idxs[1][curVec[1]] = idx;
		}
	}
}


void ofec::NBN_GridDivision::updateFitness() {


	if (!m_flag_multiThread||m_pro->hasTag(ProTag::kCSIWDN)) {
		updateFitnessThreadTask(0, m_nbn_nodes.size());
	}
	else {
		std::vector<std::thread> thrds;
		int num_task = std::thread::hardware_concurrency();
		int num_samples = m_nbn_nodes.size();
		std::vector<int> tasks;
		UTILITY::assignThreads(num_samples, num_task, tasks);
		std::pair<int, int> from_to;
		for (size_t i = 0; i < num_task; ++i) {
			from_to.first = tasks[i];
			from_to.second = tasks[i + 1];
			thrds.push_back(std::thread(
				&NBN_GridDivision::updateFitnessThreadTask, this,
				from_to.first,from_to.second));
		}
		for (auto& thrd : thrds)
			thrd.join();

	}
}


void ofec::NBN_GridDivision::updateFitnessThreadTask(int from, int to) {
	for (size_t idx(from); idx < to; ++idx) {
		m_eval_fun(*m_nbn_nodes[idx].m_sol, m_pro.get());
		m_nbn_nodes[idx].m_fitness = m_nbn_nodes[idx].m_sol->fitness();
		m_nbn_nodes[idx].m_representative = m_nbn_nodes[idx].m_sol;
	}
}

//
//void ofec::NBN_GridDivision::initNodes2D() {
//
//	if (m_pro->numVariables() == 2 && m_pro->hasTag(ofec::ProTag::kConOP)) {
//	}
//	else return;
//
//	int numSample = m_maxSample;
//	int numDim = 2;
//	int num_div = std::sqrt(numSample);
//	std::vector<double> range_div(numDim);
//	auto boundary = CAST_CONOP(m_pro)->boundary();
//	for (int idDim(0); idDim < numDim; ++idDim) {
//		auto& range_var = CAST_CONOP(m_pro)->range(idDim);
//		range_div[idDim] = (range_var.second - range_var.first) / double(num_div);
//	}
//
//	m_cur_sols.resize(numSample);
//	for (int idx(0); idx < m_cur_sols.size(); ++idx) {
//		m_cur_sols[idx].reset(new SolutionType(
//			m_problem->numObjectives(),
//			m_problem->numConstraints(),
//			m_problem->numVariables()));
//		m_cur_sols[idx]->setId(idx);
//	}
//
//	m_nbn_nodes.resize(numSample);
//	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
//		m_nbn_nodes[idx].m_sol_id = idx;
//		m_nbn_nodes[idx].m_node_idx = idx;
//	}
//
//	//for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
//	//	m_nbn_nodes[idx].m_sol = m_cur_sols[m_nbn_nodes[idx].m_sol_id];
//	//}
//	{
//		std::vector<std::vector<int>> xy_idxs(num_div);
//		for (auto& it : xy_idxs)it.resize(num_div);
//		int cur_id(0);
//
//		for (int idx(0); idx < num_div; ++idx) {
//			for (int idy(0); idy < num_div; ++idy) {
//				xy_idxs[idx][idy] = cur_id;
//
//				m_cur_sols[cur_id]->variable()[0] = boundary[0].first + idx * range_div[0];
//				m_cur_sols[cur_id]->variable()[1] = boundary[1].first + idy * range_div[1];
//				m_eval_fun(*m_cur_sols[cur_id], m_id_pro);
//
//				auto& cur_node = m_nbn_nodes[cur_id];
//				cur_node.m_vec.resize(numDim);
//				cur_node.m_vec[0] = idx;
//				cur_node.m_vec[1] = idy;
//				if (idx) {
//					int nei_id = xy_idxs[idx - 1][idy];
//					auto& nei_node = m_nbn_nodes[nei_id];
//
//					double cur_dis = m_problem->norVariableDistance(*m_cur_sols[cur_id], *m_cur_sols[nei_id]);
//					cur_node.m_neighbor.push_back(nei_id);
//					cur_node.m_neighbor_dis.push_back(cur_dis);
//					nei_node.m_neighbor.push_back(cur_id);
//					nei_node.m_neighbor_dis.push_back(cur_dis);
//
//					m_neighbors.push_back({ cur_id, nei_id });
//				}
//				if (idy) {
//					int nei_id = xy_idxs[idx][idy - 1];
//					auto& nei_node = m_nbn_nodes[nei_id];
//					double cur_dis = m_problem->norVariableDistance(*m_cur_sols[cur_id], *m_cur_sols[nei_id]);
//					cur_node.m_neighbor.push_back(nei_id);
//					cur_node.m_neighbor_dis.push_back(cur_dis);
//					nei_node.m_neighbor.push_back(cur_id);
//					nei_node.m_neighbor_dis.push_back(cur_dis);
//
//					m_neighbors.push_back({ cur_id, nei_id });
//
//				}
//				++cur_id;
//			}
//		}
//	}
//}



void ofec::NBN_GridDivision::mergeParent(Node* cur_node) {
	if (cur_node->m_parent == cur_node->m_id) return;
	std::vector<Node*> nodes;
	Node* curVisNode = cur_node;
	while (curVisNode->m_parent != curVisNode->m_id) {
		nodes.push_back(curVisNode);
		curVisNode = &m_nbn_nodes[curVisNode->m_parent];
	}
	for (auto& it : nodes) {
		it->m_parent = curVisNode->m_id;
	}
}


void ofec::NBN_GridDivision::calculate()
{

	//initNodes2D();
	//initNodesNormal();

	//std::sort(m_nbn_nodes.begin(), m_nbn_nodes.end(), [&](
	//	const node &a, const node & b
	//	) {
	//	return m_cur_sols[a.m_sol_id]->fitness() < m_cur_sols[b.m_sol_id]->fitness();
	//});

	std::vector<int> sorted_idx(m_nbn_nodes.size());
	for (int idx(0); idx < sorted_idx.size(); ++idx) {
		sorted_idx[idx] = idx;
	}
	std::sort(sorted_idx.begin(), sorted_idx.end(), [&](
		int a, int  b
		) {
		if (m_nbn_nodes[a].m_fitness == m_nbn_nodes[b].m_fitness) {
			return a < b;
		}
		else return m_nbn_nodes[a].m_fitness < m_nbn_nodes[b].m_fitness;
	});


	//std::vector<double> sorted_node_fitess(m_nbn_nodes.size());
	//for (int idx(0); idx < sorted_node_fitess.size(); ++idx) {
	//	sorted_node_fitess[idx] = m_cur_sols[idx]->fitness();
	//}


	m_cur_visited = 0;
	std::vector<int> peak_ranges;

	for (auto& sortedId : sorted_idx) {
		auto& cur_node(m_nbn_nodes[sortedId]);

		//cur_node.m_direct_parent = &cur_node;
		//double min_dis = std::numeric_limits<double>::max();
		//double cur_dis(0);
		//for (auto& nei_node : cur_node.m_neighbor) {
		//	if (nei_node->m_sol->fitness() > cur_node.m_sol->fitness()) {
		//		cur_dis = nei_node->m_sol->variableDistance(*cur_node.m_sol, m_id_pro);
		//		if (cur_dis < min_dis) {
		//			min_dis = cur_dis;
		//			cur_node.m_direct_parent = nei_node;
		//		}
		//	}
		//	//updateParent(nei_node);
		//}

		if (++m_cur_visited == 0) {
			resetVisited();
		}
		peak_ranges.clear();
		cur_node.m_visited_id = m_cur_visited;
		for (auto& nei_info : cur_node.m_neighbor_id_dis) {
			auto& nei_idx = nei_info.first;
			auto nei_node(&m_nbn_nodes[nei_idx]);
			mergeParent(nei_node);

			if (nei_node->m_parent == cur_node.m_id) {
				for (auto& son_idx : nei_node->m_better_range) {
					auto son_node(&m_nbn_nodes[son_idx]);
					mergeParent(son_node);
					if (m_nbn_nodes[son_node->m_parent].m_visited_id != m_cur_visited) {
						m_nbn_nodes[son_node->m_parent].m_visited_id = m_cur_visited;
						peak_ranges.push_back(son_node->m_parent);
					}
				}
			}
			else {
				if (m_nbn_nodes[nei_node->m_parent].m_visited_id != m_cur_visited) {
					m_nbn_nodes[nei_node->m_parent].m_visited_id = m_cur_visited;
					peak_ranges.push_back(nei_node->m_parent);
				}
			}
		}
		cur_node.m_better_range.clear();

		double min_dis = std::numeric_limits<double>::max();
		double cur_dis(0);
		double cur_dis2center(0);
		double min_dis2center = std::numeric_limits<double>::max();
		
		for (auto& nei_idx : peak_ranges) {
			auto& nei_node(m_nbn_nodes[nei_idx]);
			if (cur_node.m_sol->fitness() > nei_node.m_sol->fitness())continue;


			mergeParent(&nei_node);
			cur_dis = m_pro->norVariableDistance(*cur_node.m_sol, *nei_node.m_sol);
			cur_dis2center = m_pro->norVariableDistance(*cur_node.m_sol, *m_nbn_nodes[nei_node.m_parent].m_sol);
			
			//cur_dis = m_pro->norVariableDistance(*m_cur_sols[cur_node.m_sol_id], *m_cur_sols[m_nbn_nodes[nei_idx].m_sol_id]);
			//cur_dis = m_cur_sols[cur_node.m_sol_id]->norVariableDistance(*m_cur_sols[m_nbn_nodes[nei_idx].m_sol_id], m_id_pro);
			if (cur_dis < min_dis) {
				min_dis = cur_dis;
				min_dis2center = cur_dis2center;
				cur_node.m_parent = nei_idx;
			}			
			else if (cur_dis == min_dis) {

				if (m_random->uniform.next() < 0.5) {
					min_dis = cur_dis;
					min_dis2center = cur_dis2center;
					cur_node.m_parent = nei_idx;
				}
				//if (cur_dis2center < min_dis2center) {
				//	min_dis = cur_dis;
				//	min_dis2center = cur_dis2center;
				//	cur_node.m_parent = nei_idx;
				//}
				//else if (cur_dis2center == min_dis2center && m_random->uniform.next()<0.5) {
				//	min_dis = cur_dis;
				//	min_dis2center = cur_dis2center;
				//	cur_node.m_parent = nei_idx;
				//}
			}
		}
		cur_node.m_direct_parent = cur_node.m_parent;

		for (auto& nei_idx : peak_ranges) {
			if (nei_idx != cur_node.m_parent) {
				cur_node.m_better_range.push_back(nei_idx);
			}
		}

		{
			auto& cur_sol(cur_node.m_sol);
			auto& par_sol(m_nbn_nodes[cur_node.m_direct_parent].m_sol);
			cur_node.m_dis2parent = m_pro->norVariableDistance(*cur_sol, *par_sol);
			//cur_node.m_dis2parent = cur_sol->norVariableDistance(*par_sol, m_id_pro);
		}
	}
}



void ofec::NBN_GridDivision::updateSol(const SolBase& new_sol,int & belong_id, 
	bool flag_opt, int popIter, int popSolId, int algId ) {
	std::vector<int> cur_vec;
	int cur_id(0);
	solToVec(new_sol, cur_vec);
	vecToIdx(cur_vec, cur_id);
	auto& cur_node(m_nbn_nodes[cur_id]);
	cur_node.m_popIter = popIter;
	cur_node.m_popSolId = popSolId;
	cur_node.m_algId = algId;
	cur_node.m_flag_opt = flag_opt;
	int nearest_id(cur_id);
	double dis2nearest(m_pro->norVariableDistance(new_sol, *cur_node.m_sol));

	for (auto& nei_info : cur_node.m_neighbor_id_dis) {
		auto& nei_node(m_nbn_nodes[nei_info.first]);
		double cur_dis = m_pro->norVariableDistance(new_sol, *nei_node.m_sol);
		if (dis2nearest > cur_dis) {
			dis2nearest = cur_dis;
			nearest_id = nei_node.m_id;
		}
	}

	belong_id = nearest_id;
	std::shared_ptr<SolBase> cur_sol(m_pro->createSolution(new_sol));
	m_eval_fun(*cur_sol, m_pro.get());
	if (cur_sol->fitness() > m_nbn_nodes[nearest_id].m_fitness) {
		m_nbn_nodes[nearest_id].m_representative = cur_sol;
		m_nbn_nodes[nearest_id].m_fitness = cur_sol->fitness();
	}

	//SolutionType sol_continous(dynamic_cast<const SolutionType&>(new_sol));

	//auto& represetative_node(m_nbn_nodes[represetative_id]);
	//represetative_node.m_fitness = std::max(represetative_node.m_fitness, sol_continous.fitness());

}

//
//void ofec::NBN_GridDivision::udpate_network() {
//
//	m_belongs.resize(m_cur_sols.size());
//
//	std::vector<int> solNodeIdx(m_cur_sols.size());
//	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
//		solNodeIdx[m_nbn_nodes[idx].m_sol_id] = idx;
//	}
//	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
//		m_belongs[idx] = m_nbn_nodes[m_nbn_nodes[solNodeIdx[idx]].m_direct_parent].m_sol_id;
//	}
//
//	m_dis2parent.resize(m_cur_sols.size());
//	for (int idx(0); idx < m_cur_sols.size(); ++idx) {
//		if (m_belongs[idx] == idx) {
//			m_dis2parent[idx] = std::numeric_limits<double>::max();
//		}
//		else {
//			m_dis2parent[idx] = m_cur_sols[idx]->variableDistance(*m_cur_sols[m_belongs[idx]], m_id_pro);
//		}
//	}
//}

void ofec::NBN_GridDivision::getNearestBetterNetworkShareMemory(
	std::vector<std::shared_ptr<SolBase>>& sols,
	std::vector<double>& fitness,
	std::vector<int>& belong, std::vector<bool>& flagOpt)
{
	//udpate_network();
	//sols = m_cur_sols;
	sols.resize(m_nbn_nodes.size());
	fitness.resize(m_nbn_nodes.size());
	belong.resize(m_nbn_nodes.size());
	flagOpt.resize(m_nbn_nodes.size());

	for (int idx(0); idx < sols.size(); ++idx) {
		sols[idx] = m_nbn_nodes[idx].m_sol;
		//sols[idx].reset(new SolutionType(*m_nbn_nodes[idx].m_sol));
		fitness[idx] = m_nbn_nodes[idx].m_fitness;
		belong[idx] = m_nbn_nodes[idx].m_direct_parent;
		flagOpt[idx] = m_nbn_nodes[idx].m_flag_opt;
	}


	//for (int idx(0); idx < sols.size(); ++idx) {
	//	sols[idx].reset(new SolutionType(*m_nbn_nodes[idx].m_sol));
	//}
	//belong = m_belongs;
	//fitness.resize(sols.size());
	//for (int idx(0); idx < sols.size(); ++idx) {
	//	fitness[idx] = sols[idx]->fitness();
	//}
}

void  ofec::NBN_GridDivision::getNearestBetterNetworkShareMemory(
	std::vector<std::shared_ptr<SolBase>>& sols,
	std::vector<std::shared_ptr<SolBase>>& representative,
	std::vector<double>& fitness,
	std::vector<int>& belong,
	std::vector<double>& dis2par
) {
	sols.resize(m_nbn_nodes.size());
	representative.resize(m_nbn_nodes.size());
	fitness.resize(m_nbn_nodes.size());
	belong.resize(m_nbn_nodes.size());
	dis2par.resize(m_nbn_nodes.size());
	for (int idx(0); idx < sols.size(); ++idx) {
		sols[idx]= m_nbn_nodes[idx].m_sol;
		representative[idx] = m_nbn_nodes[idx].m_representative;
		fitness[idx] = m_nbn_nodes[idx].m_fitness;
		belong[idx] = m_nbn_nodes[idx].m_direct_parent;
		dis2par[idx] = m_nbn_nodes[idx].m_dis2parent;
	}
}

void ofec::NBN_GridDivision::getNearestBetterNetworkShareMemory(
	std::vector<std::shared_ptr<SolBase>>& sols,
	std::vector<std::shared_ptr<SolBase>>& representative,
	std::vector<double>& fitness,
	std::vector<int>& belong,
	std::vector<double>& dis2par,
	std::vector<int> popIters,
	std::vector<int> popSolIds,
	std::vector<int> algIds
) {
	sols.resize(m_nbn_nodes.size());
	representative.resize(m_nbn_nodes.size());
	fitness.resize(m_nbn_nodes.size());
	belong.resize(m_nbn_nodes.size());
	dis2par.resize(m_nbn_nodes.size());
	popIters.resize(m_nbn_nodes.size());
	popSolIds.resize(m_nbn_nodes.size());
	algIds.resize(m_nbn_nodes.size());
	for (int idx(0); idx < sols.size(); ++idx) {
		sols[idx] = m_nbn_nodes[idx].m_sol;
		representative[idx] = m_nbn_nodes[idx].m_representative;
		fitness[idx] = m_nbn_nodes[idx].m_fitness;
		belong[idx] = m_nbn_nodes[idx].m_direct_parent;
		dis2par[idx] = m_nbn_nodes[idx].m_dis2parent;
		popIters[idx] = m_nbn_nodes[idx].m_popIter;
		popSolIds[idx] = m_nbn_nodes[idx].m_popSolId;
		algIds[idx] = m_nbn_nodes[idx].m_algId;
		
	}
}



void ofec::NBN_GridDivision::getSharedNBN(
	std::vector<std::shared_ptr<SolBase>>& sols,
	std::vector<double>& fitness,
	std::vector<int>& belong,
	std::vector<double>& dis2parent, std::vector<bool>& flagOpt
	)const{
	

	sols.resize(m_nbn_nodes.size());
	fitness.resize(m_nbn_nodes.size());
	belong.resize(m_nbn_nodes.size());
	dis2parent.resize(m_nbn_nodes.size());
	flagOpt.resize(m_nbn_nodes.size());
	for (int idx(0); idx < sols.size(); ++idx) {
		sols[idx] = m_nbn_nodes[idx].m_representative;
		fitness[idx] = m_nbn_nodes[idx].m_fitness;
		belong[idx] = m_nbn_nodes[idx].m_direct_parent;
		dis2parent[idx] = m_nbn_nodes[idx].m_dis2parent;
		flagOpt[idx] = m_nbn_nodes[idx].m_flag_opt;
	}
}


void ofec::NBN_GridDivision::getSharedNBN(
	std::vector<std::shared_ptr<SolBase>>& sols,
	std::vector<double>& fitness,
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	std::vector<int> popIters,
	std::vector<int> popSolIds,
	std::vector<int> algIds
)const {
	sols.resize(m_nbn_nodes.size());
	fitness.resize(m_nbn_nodes.size());
	belong.resize(m_nbn_nodes.size());
	dis2parent.resize(m_nbn_nodes.size());

	popIters.resize(m_nbn_nodes.size());
	popSolIds.resize(m_nbn_nodes.size());
	algIds.resize(m_nbn_nodes.size());
	for (int idx(0); idx < sols.size(); ++idx) {
		sols[idx] = m_nbn_nodes[idx].m_representative;
		fitness[idx] = m_nbn_nodes[idx].m_fitness;
		belong[idx] = m_nbn_nodes[idx].m_direct_parent;
		dis2parent[idx] = m_nbn_nodes[idx].m_dis2parent;

		popIters[idx] = m_nbn_nodes[idx].m_popIter;
		popSolIds[idx] = m_nbn_nodes[idx].m_popSolId;
		algIds[idx] = m_nbn_nodes[idx].m_algId;
	}
}



//
//void ofec::NBN_grid_sample::outputData(SampleNearestBetterNetworkRecord& record) {
//	record.m_samples.resize(m_cur_sols.size());
//	record.m_sorted_ids.resize(m_cur_sols.size());
//	record.m_parent_ids.resize(m_cur_sols.size());
//	record.m_neighbor_ids.resize(m_cur_sols.size());
//	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
//		m_nbn_nodes[idx].m_node_idx = idx;
//	}
//	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
//		record.m_samples[idx] = m_cur_sols[m_nbn_nodes[idx].m_sol_id];
//		record.m_sorted_ids[idx] = idx;
//		record.m_parent_ids[idx] = m_nbn_nodes[m_nbn_nodes[idx].m_direct_parent].m_node_idx;
//		for (auto& nei_idx : m_nbn_nodes[idx].m_neighbor) {
//			record.m_neighbor_ids[idx].push_back(nei_idx);
//		}
//	}
//}

//
//void ofec::NBN_GridSample::getTree() {
//	m_tree.resize(m_nbn_nodes.size());
//	int tree_head(0);
//	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
//		m_nbn_nodes[idx].m_node_idx = idx;
//		if (m_nbn_nodes[idx].m_direct_parent == m_nbn_nodes[idx].m_node_idx) {
//			m_head = &m_tree[idx];
//			++tree_head;
//		}
//	}
//	for (int idx(0); idx < m_tree.size(); ++idx) {
//		m_tree[idx].m_sol = m_cur_sols[idx];
//
//		m_tree[m_nbn_nodes[idx].m_direct_parent].m_leaves.push_back(&m_tree[idx]);
//	}
//
//}
//





void ofec::NBN_GridDivision::getNearestBetterNetworkShareMemoryFilter(
	std::vector<std::shared_ptr<SolBase>>& sols,
	std::vector<double>& fitness,
	std::vector<int>& belong) {

	double filterDis = m_filterDis;
	std::vector<bool> flagOpt;
	getNearestBetterNetworkShareMemory(sols, fitness, belong, flagOpt);
	std::vector<std::vector<int>> sons(belong.size());

	for (int idx(0); idx < belong.size(); ++idx) {
		if (belong[idx] != idx) {
			sons[belong[idx]].push_back(idx);
		}
	}

	std::vector<bool> flag_del(belong.size(), false);
	std::vector<int> sortedIds(sols.size());
	for (int idx(0); idx < sortedIds.size(); ++idx) {
		sortedIds[idx] = idx;
	}

	std::sort(sortedIds.begin(), sortedIds.end(), [&](int a, int b) {
		return sols[a]->fitness() < sols[b]->fitness();
	});

	for (auto& curId : sortedIds) {
		if (sons[curId].size() != 0 && curId != belong[curId]) {
			auto& curSol(sols[curId]);
			auto& parSol(sols[belong[curId]]);
			double maxDis(0);
			for (auto& sonId : sons[curId]) {
				maxDis = std::max(maxDis, m_pro->norVariableDistance(*parSol, *sols[sonId]));
				//parSol->norVariableDistance(*sols[sonId], m_id_pro));
//				maxDis = std::max(maxDis, parSol->norVariableDistance(*sols[sonId], m_id_pro));
			}
			if (maxDis < filterDis) {
				flag_del[curId] = true;
				auto& par_son(sons[belong[curId]]);
				for (int sonIter(0); sonIter < sons[curId].size(); ++sonIter) {
					int sonId = sons[curId][sonIter];
					belong[sonId] = belong[curId];
					par_son.push_back(sonId);
				}
				//for (auto& sonId : sons[curId]) {

				//}
			}
		}
	}

	std::vector<std::shared_ptr<SolBase>> filter_sols;
	std::vector<int> newToOld;
	std::vector<int> oldToNew(sols.size(), -1);
	std::vector<int> newBelong;

	{
		int newIdx(0);
		for (int idx(0); idx < sols.size(); ++idx) {
			if (!flag_del[idx]) {
				newToOld.push_back(idx);
				oldToNew[idx] = newIdx++;
				filter_sols.push_back(sols[idx]);
			}
		}

		newBelong.resize(newIdx);
		for (int nIdx(0); nIdx < newBelong.size(); ++nIdx) {
			newBelong[nIdx] = oldToNew[belong[newToOld[nIdx]]];
		}

	}

	sols = filter_sols;
	belong = newBelong;
	fitness.resize(sols.size());
	for (int idx(0); idx < fitness.size(); ++idx) {
		fitness[idx] = sols[idx]->fitness();
	}


}

//
//
//int ofec::NBN_GridDivision::addRandomSol(const std::shared_ptr<SolutionType>& new_sol, int to_idx) {
//	m_cur_sols.push_back(new_sol);
//	int cur_idx(m_nbn_nodes.size());
//	m_nbn_nodes.push_back(Node());
//	auto& cur_node(m_nbn_nodes.back());
//	cur_node.clear();
//	cur_node.init(cur_idx);
//	//cur_node.m_node_idx = cur_node.m_sol_id = cur_idx;
//	auto& cur_sol(m_cur_sols.back());
//
//	double cur_dis(0), before_dis(0);
//	for (int nei_idx(0); nei_idx < to_idx; ++nei_idx) {
//		auto& nei_node(m_nbn_nodes[nei_idx]);
//		auto& nei_sol(m_cur_sols[nei_node.m_sol_id]);
//
//		cur_dis = m_pro->norVariableDistance(*cur_sol, *nei_sol);
//
//		//cur_dis = cur_sol->norVariableDistance(*nei_sol, m_id_pro);
//		if (nei_sol->fitness() > cur_sol->fitness()) {
//			if (cur_dis < cur_node.m_dis2parent) {
//				cur_node.m_dis2parent = cur_dis;
//				cur_node.m_direct_parent = nei_idx;
//			}
//		}
//		else if (nei_sol->fitness() < cur_sol->fitness()) {
//			if (cur_dis < nei_node.m_dis2parent) {
//				nei_node.m_dis2parent = cur_dis;
//				nei_node.m_direct_parent = cur_idx;
//			}
//		}
//	}
//
//	return cur_idx;
//
//}
//
//
//int ofec::NBN_GridSample::addRandomSol(const std::shared_ptr<SolutionType>& new_sol) {
//	m_cur_sols.push_back(new_sol);
//	int cur_idx(m_nbn_nodes.size());
//	m_nbn_nodes.push_back(Node());
//	auto& cur_node(m_nbn_nodes.back());
//	auto& cur_sol(m_cur_sols.back());
//	cur_node.clear();
//	cur_node.init(cur_idx);
//	//cur_node.m_node_idx = cur_node.m_sol_id = cur_idx;
//
//
//	double cur_dis(0), before_dis(0);
//	for (int nei_idx(0); nei_idx < cur_idx; ++nei_idx) {
//		auto& nei_node(m_nbn_nodes[nei_idx]);
//		auto& nei_sol(m_cur_sols[nei_node.m_sol_id]);
//		cur_dis = m_pro->norVariableDistance(*cur_sol, *nei_sol);
//		//		cur_dis = cur_sol->norVariableDistance(*nei_sol, m_id_pro);
//		if (nei_sol->fitness() > cur_sol->fitness()) {
//			if (cur_dis < cur_node.m_dis2parent) {
//				cur_node.m_dis2parent = cur_dis;
//				cur_node.m_direct_parent = nei_idx;
//			}
//		}
//		else if (nei_sol->fitness() < cur_sol->fitness()) {
//			if (cur_dis < nei_node.m_dis2parent) {
//				nei_node.m_dis2parent = cur_dis;
//				nei_node.m_direct_parent = cur_idx;
//			}
//		}
//	}
//
//	return cur_idx;
//
//}
//
//void ofec::NBN_GridSample::insertOpt() {
//	std::shared_ptr<SolutionType> cur_sol;
//	int num_vars = m_pro->numVariables();
//	int num_objs = m_pro->numObjectives();
//	int num_cons = m_pro->numConstraints();
//
//	auto& optBase(m_pro->optBase());
//
//	m_marker_idxs.clear();
//	int to_idx(m_nbn_nodes.size());
//	int node_idx(0);
//	//std::vector<std::shared_ptr<SolutionType>> opt_sols;
//	for (size_t idx(0); idx < optBase.numberVariables(); ++idx) {
//		cur_sol.reset(new SolutionType(num_objs, num_cons, num_vars));
//		cur_sol->variable() = dynamic_cast<const SolutionType&>(optBase.variable(idx)).variable();
//		m_eval_fun(*cur_sol, m_id_pro);
//		m_marker_idxs.push_back(addRandomSol(cur_sol, to_idx));
//		//opt_sols.push_back(cur_sol);
//		//m_nbn_network.addRandomSol(cur_sol);
//	}
//
//	for (auto& cur_idx : m_marker_idxs) {
//		auto& cur_node(m_nbn_nodes[cur_idx]);
//		auto& cur_sol(m_cur_sols[cur_node.m_sol_id]);
//		for (auto& nei_idx : m_marker_idxs) {
//			if (cur_idx != nei_idx) {
//				auto& nei_node(m_nbn_nodes[nei_idx]);
//				auto& nei_sol(m_cur_sols[nei_node.m_sol_id]);
//				cur_node.m_dis2parent = std::min(cur_node.m_dis2parent,
//					m_pro->norVariableDistance(*cur_sol, *nei_sol));
//				//cur_sol->norVariableDistance(*nei_sol, m_id_pro));
//
////				cur_node.m_dis2parent = std::min(cur_node.m_dis2parent, cur_sol->norVariableDistance(*nei_sol, m_id_pro));
//			}
//		}
//	}
//
//
//	std::cout << "dynamic sample alg opt numbers\t" << optBase.numberVariables() << std::endl;
//	//#ifdef  OFEC_DEMO
//	//
//	//	auto& markderIdxs(ofec_demo::bufferAlgNearestBetterNetwork::ms_markerPonits);
//	//	markderIdxs.clear();
//	//	for (size_t idx(0); idx < optBase.numberVariables(); ++idx) {
//	//		markderIdxs.push_back(idx);
//	//	}
//	//
//	//#endif //  OFEC_DEMO
//}
