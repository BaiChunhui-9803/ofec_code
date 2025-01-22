#include "sample2d_graph_algorithm.h"
#include "../../../core/problem/continuous/continuous.h"

void ofec::Sample2D_Graph_Algorithm::init(Problem *pro, Random *rnd, const std::function<void(Solution<>& sol, Problem *pro)>& eval_fun)
{
	m_random.get() = rnd;
	m_problem.get() = pro;
	m_eval_fun = eval_fun;
	m_dim = m_problem->numberVariables();
}


void ofec::Sample2D_Graph_Algorithm::idxToVec(int idx, std::vector<int>& cur) {
	
	for (int idDim(m_dim - 1); idDim >= 0; --idDim) {
		cur[idDim] = idx % m_dim_div[idDim];
		idx /= m_dim_div[idDim];
	}
	
}
void ofec::Sample2D_Graph_Algorithm::vecToIdx(const std::vector<int>& cur, int& idx) {
	//for(int idDim)
	idx = 0;
	for (int idDim(0); idDim < m_dim; ++idDim) {
		idx *= m_dim_div[idDim];
		idx += cur[idDim];
	}

	//for (int idDim(m_dim - 1); idDim >= 0; --idDim) {
	////	cur[idDim] = idx % m_dim_div[idDim];
	////	idx /= m_dim_div[idDim];
	//}
}
bool ofec::Sample2D_Graph_Algorithm::judgeFeasible(const std::vector<int>& cur) {
	for (int idx(0); idx < m_dim; ++idx) {
		if (cur[idx] >= 0 && cur[idx] < m_dim_div[idx]) {

		}
		else return false;
	}
	return true;
}

void ofec::Sample2D_Graph_Algorithm::initNodesNormal() {
	
	//std::vector<Node*> nodes;
	
	int numSample_int = m_maxSample;
	m_dim_div.resize(m_dim);
	int dim_div = exp(log(double(numSample_int))/ double(m_dim));
	numSample_int = 1;
	for (auto& it : m_dim_div) {
		it = dim_div;
		numSample_int *= it;
	}

	m_cur_sols.resize(numSample_int);
	for (int idx(0); idx < m_cur_sols.size(); ++idx) {
		m_cur_sols[idx].reset(new SolutionType(
			m_problem->numberObjectives(),
			m_problem->numberConstraints(),
			m_problem->numberVariables()));
		m_cur_sols[idx]->setId(idx);
	}
	m_nbn_nodes.resize(numSample_int);

	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
		m_nbn_nodes[idx].m_sol_id = idx;
		m_nbn_nodes[idx].m_node_idx = idx;
	}
	//for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
	//	m_nbn_nodes[idx].m_sol = m_cur_sols[m_nbn_nodes[idx].m_sol_id];
	//}

	std::vector<double> range_div(m_dim);
	auto boundary = CAST_CONOP(m_problem.get())->boundary();
	for (int idDim(0); idDim < m_dim; ++idDim) {
		auto& range_var = CAST_CONOP(m_problem.get())->range(idDim);
		range_div[idDim] = (range_var.second - range_var.first) / double(m_dim_div[idDim]);
	}

	std::vector<int> curVec(m_dim,0);
	int curIdx(0);

	// set solutions
	for (size_t idx(0); idx < numSample_int; ++idx) {
		idxToVec(idx, curVec);
		vecToIdx(curVec, curIdx);
		m_nbn_nodes[idx].m_vec = curVec;
		if (idx != curIdx) {
			THROW("error at wrong vec to idx \n");
			//std::cout << "error" << endl;
		}
		for (int idDim(0); idDim < m_dim; ++idDim) {
			m_cur_sols[idx]->variable()[idDim] = boundary[idDim].first + curVec[idDim] * range_div[idDim];
		}
		m_eval_fun(*m_cur_sols[idx], m_problem.get());
	}

	// set neighbor sols
	std::vector<int> neiVec(m_dim, 0);
	int neiIdx(0);
	for (int idx(0); idx < numSample_int; ++idx) {
		auto& cur_node = m_nbn_nodes[idx];
		idxToVec(idx, curVec);
		vecToIdx(curVec, curIdx);
		neiVec = curVec;
		for (int idDim(0); idDim < m_dim; ++idDim) {
			if (--neiVec[idDim] >= 0) {
				vecToIdx(neiVec, neiIdx);
				auto& nei_node = m_nbn_nodes[neiIdx];

				double cur_dis = m_problem->normalizedVariableDistance(*m_cur_sols[idx], *m_cur_sols[neiIdx]);
				cur_node.m_neighbor.push_back(neiIdx);
				cur_node.m_neighbor_dis.push_back(cur_dis);
				nei_node.m_neighbor.push_back(curIdx);
				nei_node.m_neighbor_dis.push_back(cur_dis);

				m_neighbors.push_back({ idx, neiIdx });
			}
			++neiVec[idDim];
		}
	}

	


}
void ofec::Sample2D_Graph_Algorithm::initNodes2D() {

	if (m_problem->numberVariables() == 2 && m_problem->hasTag(ofec::ProblemTag::kConOP)) {
	}
	else return;

	int numSample = m_maxSample;
	int numDim = 2;
	int num_div = std::sqrt(numSample);
	std::vector<double> range_div(numDim);
	auto boundary = CAST_CONOP(m_problem.get())->boundary();
	for (int idDim(0); idDim < numDim; ++idDim) {
		auto& range_var = CAST_CONOP(m_problem.get())->range(idDim);
		range_div[idDim] = (range_var.second - range_var.first) / double(num_div);
	}

	m_cur_sols.resize(numSample);
	for (int idx(0); idx < m_cur_sols.size(); ++idx) {
		m_cur_sols[idx].reset(new SolutionType(
			m_problem->numberObjectives(),
			m_problem->numberConstraints(),
			m_problem->numberVariables()));
		m_cur_sols[idx]->setId(idx);
	}

	m_nbn_nodes.resize(numSample);
	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
		m_nbn_nodes[idx].m_sol_id = idx;
		m_nbn_nodes[idx].m_node_idx = idx;
	}

	//for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
	//	m_nbn_nodes[idx].m_sol = m_cur_sols[m_nbn_nodes[idx].m_sol_id];
	//}
	{
		std::vector<std::vector<int>> xy_idxs(num_div);
		for (auto& it : xy_idxs)it.resize(num_div);
		int cur_id(0);

		for (int idx(0); idx < num_div; ++idx) {
			for (int idy(0); idy < num_div; ++idy) {
				xy_idxs[idx][idy] = cur_id;

				m_cur_sols[cur_id]->variable()[0] = boundary[0].first + idx * range_div[0];
				m_cur_sols[cur_id]->variable()[1] = boundary[1].first + idy * range_div[1];
				m_eval_fun(*m_cur_sols[cur_id], m_problem.get());

				auto& cur_node = m_nbn_nodes[cur_id];
				cur_node.m_vec.resize(numDim);
				cur_node.m_vec[0] = idx;
				cur_node.m_vec[1] = idy;
				if (idx) {
					int nei_id = xy_idxs[idx - 1][idy];
					auto& nei_node = m_nbn_nodes[nei_id];

					double cur_dis = m_problem->normalizedVariableDistance(*m_cur_sols[cur_id], *m_cur_sols[nei_id]);
					cur_node.m_neighbor.push_back(nei_id);
					cur_node.m_neighbor_dis.push_back(cur_dis);
					nei_node.m_neighbor.push_back(cur_id);
					nei_node.m_neighbor_dis.push_back(cur_dis);

					m_neighbors.push_back({ cur_id, nei_id });
				}
				if (idy) {
					int nei_id = xy_idxs[idx][idy - 1];
					auto& nei_node = m_nbn_nodes[nei_id];
					double cur_dis = m_problem->normalizedVariableDistance(*m_cur_sols[cur_id], *m_cur_sols[nei_id]);
					cur_node.m_neighbor.push_back(nei_id);
					cur_node.m_neighbor_dis.push_back(cur_dis);
					nei_node.m_neighbor.push_back(cur_id);
					nei_node.m_neighbor_dis.push_back(cur_dis);

					m_neighbors.push_back({ cur_id, nei_id });

				}
				++cur_id;
			}
		}
	}
}

void ofec::Sample2D_Graph_Algorithm::setParents() {
	for (auto& it : m_nbn_nodes) it.m_parent = it.m_node_idx;
	//for(int idx(0);idx<)
}

void ofec::Sample2D_Graph_Algorithm::updateParent(node * cur_node) {
	if (cur_node->m_parent == cur_node->m_node_idx) return;	
	std::vector<node*> nodes;
	node* curVisNode = cur_node;
	while (curVisNode->m_parent != curVisNode->m_node_idx) {
		nodes.push_back(curVisNode);
		curVisNode = &m_nbn_nodes[curVisNode->m_parent];
	}
	for (auto& it : nodes) {
		it->m_parent = curVisNode->m_node_idx;
	}
}


void ofec::Sample2D_Graph_Algorithm::calculate()
{

	//initNodes2D();
	initNodesNormal();

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
		int a ,int  b
		) {
		return m_cur_sols[a]->fitness() < m_cur_sols[b]->fitness();
	});


	setParents();
	std::vector<double> sorted_node_fitess(m_cur_sols.size());
	for (int idx(0); idx < sorted_node_fitess.size(); ++idx) {
		sorted_node_fitess[idx] =m_cur_sols[idx]->fitness();
	}


	m_cur_visited = 0;
	std::vector<int> peak_ranges;

	for (auto& sortedId : sorted_idx) {
		auto& cur_node(m_nbn_nodes[sortedId]);

		//cur_node.m_direct_parent = &cur_node;
		//double min_dis = std::numeric_limits<double>::max();
		//double cur_dis(0);
		//for (auto& nei_node : cur_node.m_neighbor) {
		//	if (nei_node->m_sol->fitness() > cur_node.m_sol->fitness()) {
		//		cur_dis = nei_node->m_sol->variableDistance(*cur_node.m_sol, m_problem.get());
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
		for (auto& nei_idx : cur_node.m_neighbor) {
			auto nei_node(&m_nbn_nodes[nei_idx]);
			updateParent(nei_node);

			if (nei_node->m_parent == cur_node.m_node_idx) {
				for (auto& son_idx : nei_node->m_better_range) {
					auto son_node(&m_nbn_nodes[son_idx]);
					updateParent(son_node);
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
		for (auto& nei_idx : peak_ranges) {
			cur_dis = m_problem->normalizedVariableDistance(*m_cur_sols[cur_node.m_sol_id], *m_cur_sols[m_nbn_nodes[nei_idx].m_sol_id]);
			//cur_dis = m_cur_sols[cur_node.m_sol_id]->normalizedVariableDistance(*m_cur_sols[m_nbn_nodes[nei_idx].m_sol_id], m_problem.get());
			if (cur_dis < min_dis) {
				min_dis = cur_dis;
				if (cur_node.m_parent != cur_node.m_node_idx) {
					cur_node.m_better_range.push_back(cur_node.m_parent);
				}
				cur_node.m_parent = nei_idx;
			}
			else if (cur_dis == min_dis && m_random->uniform.next() < 0.5) {
				min_dis = cur_dis;
				if (cur_node.m_parent != cur_node.m_node_idx) {
					cur_node.m_better_range.push_back(cur_node.m_parent);
				}
				cur_node.m_parent = nei_idx;
			}
			else {
				cur_node.m_better_range.push_back(nei_idx);
			}
		}
		cur_node.m_direct_parent = cur_node.m_parent;

		{
			auto& cur_sol(m_cur_sols[cur_node.m_sol_id]);
			auto& par_sol(m_cur_sols[m_nbn_nodes[cur_node.m_direct_parent].m_sol_id]);
			cur_node.m_dis2parent = m_problem->normalizedVariableDistance(*cur_sol, *par_sol);
			//cur_node.m_dis2parent = cur_sol->normalizedVariableDistance(*par_sol, m_problem.get());
		}
	}
}


void ofec::Sample2D_Graph_Algorithm::udpate_network() {
	
	m_belongs.resize(m_cur_sols.size());
	std::vector<int> solNodeIdx(m_cur_sols.size());
	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
		solNodeIdx[m_nbn_nodes[idx].m_sol_id] = idx;
	}
	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
		m_belongs[idx] = m_nbn_nodes[m_nbn_nodes[solNodeIdx[idx]].m_direct_parent].m_sol_id;
	}
}

void ofec::Sample2D_Graph_Algorithm::getNearestBetterNetworkShareMemory(
	std::vector<std::shared_ptr<SolutionBase>>& sols,
	std::vector<double>& fitness, 
	std::vector<int>& belong)
{
	udpate_network();
	//sols = m_cur_sols;
	sols.resize(m_cur_sols.size());
	for (int idx(0); idx < sols.size(); ++idx) {
		sols[idx].reset(new SolutionType(*m_cur_sols[idx]));
	}
	belong = m_belongs;
	fitness.resize(sols.size());
	for (int idx(0); idx < sols.size(); ++idx) {
		fitness[idx] = sols[idx]->fitness();
	}
}

void  ofec::Sample2D_Graph_Algorithm::getNearestBetterNetworkShareMemory(
	std::vector<std::shared_ptr<SolutionBase>>& sols,
	std::vector<double>& fitness,
	std::vector<int>& belong,
	std::vector<double>& dis2par
) {

}





void ofec::Sample2D_Graph_Algorithm::outputData(SampleNearestBetterNetworkRecord& record) {
	record.m_samples.resize(m_cur_sols.size());
	record.m_sorted_ids.resize(m_cur_sols.size());
	record.m_parent_ids.resize(m_cur_sols.size());
	record.m_neighbor_ids.resize(m_cur_sols.size());
	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
		m_nbn_nodes[idx].m_node_idx = idx;	
	}
	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
		record.m_samples[idx] = m_cur_sols[m_nbn_nodes[idx].m_sol_id];
		record.m_sorted_ids[idx] = idx;
		record.m_parent_ids[idx] = m_nbn_nodes[m_nbn_nodes[idx].m_direct_parent].m_node_idx;
		for (auto& nei_idx : m_nbn_nodes[idx].m_neighbor) {
			record.m_neighbor_ids[idx].push_back(nei_idx);
		}
	}
}


void ofec::Sample2D_Graph_Algorithm::getTree() {
	m_tree.resize(m_nbn_nodes.size());
	int tree_head(0);
	for (int idx(0); idx < m_nbn_nodes.size(); ++idx) {
		m_nbn_nodes[idx].m_node_idx = idx;
		if (m_nbn_nodes[idx].m_direct_parent == m_nbn_nodes[idx].m_node_idx) {
			m_head = &m_tree[idx];
			++tree_head;
		}
	}
	for (int idx(0); idx < m_tree.size(); ++idx) {
		m_tree[idx].m_sol = m_cur_sols[idx];
		
		m_tree[m_nbn_nodes[idx].m_direct_parent].m_leaves.push_back(&m_tree[idx]);
	}
	
}






void ofec::Sample2D_Graph_Algorithm::getNearestBetterNetworkShareMemoryFilter(
	std::vector<std::shared_ptr<SolutionBase>>& sols,
	std::vector<double>& fitness,
	std::vector<int>& belong) {
	
	double filterDis = m_filterDis;
	udpate_network();
	getNearestBetterNetworkShareMemory(sols, fitness, belong);
	std::vector<std::vector<int>> sons(belong.size());

	for (int idx(0); idx < belong.size(); ++idx) {
		if (belong[idx] != idx) {
			sons[belong[idx]].push_back(idx);
		}
	}

	std::vector<bool> flag_del(belong.size(),false);
	std::vector<int> sortedIds(sols.size());
	for (int idx(0); idx < sortedIds.size(); ++idx) {
		sortedIds[idx] = idx;
	}

	std::sort(sortedIds.begin(), sortedIds.end(), [&](int a,int b) {
		return sols[a]->fitness() < sols[b]->fitness();
	});

	for (auto& curId : sortedIds) {
		if (sons[curId].size() != 0&&curId!=belong[curId]) {
			auto& curSol(sols[curId]);
			auto& parSol(sols[belong[curId]]);
			double maxDis(0);
			for (auto& sonId : sons[curId]) {
				maxDis = std::max(maxDis, m_problem->normalizedVariableDistance(*parSol, *sols[sonId]));
					//parSol->normalizedVariableDistance(*sols[sonId], m_problem.get()));
//				maxDis = std::max(maxDis, parSol->normalizedVariableDistance(*sols[sonId], m_problem.get()));
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
	
	std::vector<std::shared_ptr<SolutionBase>> filter_sols;
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



int ofec::Sample2D_Graph_Algorithm::addRandomSol(const std::shared_ptr<SolutionType>& new_sol, int to_idx) {
	m_cur_sols.push_back(new_sol);
	int cur_idx(m_nbn_nodes.size());
	m_nbn_nodes.push_back(node());
	auto& cur_node(m_nbn_nodes.back());
	cur_node.clear();
	cur_node.init(cur_idx);
	//cur_node.m_node_idx = cur_node.m_sol_id = cur_idx;
	auto& cur_sol(m_cur_sols.back());

	double cur_dis(0), before_dis(0);
	for (int nei_idx(0); nei_idx < to_idx; ++nei_idx) {
		auto& nei_node(m_nbn_nodes[nei_idx]);
		auto& nei_sol(m_cur_sols[nei_node.m_sol_id]);

		cur_dis = m_problem->normalizedVariableDistance(*cur_sol, *nei_sol);

		//cur_dis = cur_sol->normalizedVariableDistance(*nei_sol, m_problem.get());
		if (nei_sol->fitness() > cur_sol->fitness()) {
			if (cur_dis < cur_node.m_dis2parent) {
				cur_node.m_dis2parent = cur_dis;
				cur_node.m_direct_parent = nei_idx;
			}
		}
		else if (nei_sol->fitness() < cur_sol->fitness()) {
			if (cur_dis < nei_node.m_dis2parent) {
				nei_node.m_dis2parent = cur_dis;
				nei_node.m_direct_parent = cur_idx;
			}
		}
	}

	return cur_idx;

}


int ofec::Sample2D_Graph_Algorithm::addRandomSol(const std::shared_ptr<SolutionType>& new_sol) {
	m_cur_sols.push_back(new_sol);
	int cur_idx(m_nbn_nodes.size());
	m_nbn_nodes.push_back(node());
	auto& cur_node(m_nbn_nodes.back());
	auto& cur_sol(m_cur_sols.back());
	cur_node.clear();
	cur_node.init(cur_idx);
	//cur_node.m_node_idx = cur_node.m_sol_id = cur_idx;


	double cur_dis(0), before_dis(0);
	for (int nei_idx(0); nei_idx < cur_idx; ++nei_idx) {
		auto& nei_node(m_nbn_nodes[nei_idx]);
		auto& nei_sol(m_cur_sols[nei_node.m_sol_id]);
		cur_dis = m_problem->normalizedVariableDistance(*cur_sol,*nei_sol);
//		cur_dis = cur_sol->normalizedVariableDistance(*nei_sol, m_problem.get());
		if (nei_sol->fitness() > cur_sol->fitness()) {
			if (cur_dis < cur_node.m_dis2parent) {
				cur_node.m_dis2parent = cur_dis;
				cur_node.m_direct_parent = nei_idx;
			}
		}
		else if (nei_sol->fitness() < cur_sol->fitness()) {
			if (cur_dis < nei_node.m_dis2parent) {
				nei_node.m_dis2parent = cur_dis;
				nei_node.m_direct_parent = cur_idx;
			}
		}
	}

	return cur_idx;

}

void ofec::Sample2D_Graph_Algorithm::insertOpt() {
	std::shared_ptr<SolutionType> cur_sol;
	int num_vars = m_problem->numberVariables();
	int number_objectives = m_problem->numberObjectives();
	int num_cons = m_problem->numberConstraints();

	auto& optBase(m_problem->optBase());

	m_marker_idxs.clear();
	int to_idx(m_nbn_nodes.size());
	int node_idx(0);
	//std::vector<std::shared_ptr<SolutionType>> opt_sols;
	for (size_t idx(0); idx < optBase.numberVariables(); ++idx) {
		cur_sol.reset(new SolutionType(number_objectives, num_cons, num_vars));
		cur_sol->variable() = dynamic_cast<const SolutionType&>(optBase.variable(idx)).variable();
		m_eval_fun(*cur_sol, m_problem.get());
		m_marker_idxs.push_back(addRandomSol(cur_sol, to_idx));
		//opt_sols.push_back(cur_sol);
		//m_nbn_network.addRandomSol(cur_sol);
	}

	for (auto& cur_idx : m_marker_idxs) {
		auto& cur_node(m_nbn_nodes[cur_idx]);
		auto& cur_sol(m_cur_sols[cur_node.m_sol_id]);
		for (auto& nei_idx : m_marker_idxs) {
			if (cur_idx != nei_idx) {
				auto& nei_node(m_nbn_nodes[nei_idx]);
				auto& nei_sol(m_cur_sols[nei_node.m_sol_id]);
				cur_node.m_dis2parent = std::min(cur_node.m_dis2parent,
					m_problem->normalizedVariableDistance(*cur_sol, *nei_sol));
					//cur_sol->normalizedVariableDistance(*nei_sol, m_problem.get()));

//				cur_node.m_dis2parent = std::min(cur_node.m_dis2parent, cur_sol->normalizedVariableDistance(*nei_sol, m_problem.get()));
			}
		}
	}


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
