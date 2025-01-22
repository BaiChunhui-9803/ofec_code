#include "nbn_alg_com_gl_pop.h"

int ofec::PopGL_NBN_COM_ALG::evolve(Problem* pro, Algorithm* alg, Random* rnd)
{


	auto tag = evaluate(pro, alg);
	int numCity = pro->numVariables();
	std::vector<int> visited(numCity);
//	std::vector<double> pro;
	int from = 0;
	std::function< double(const int& cur_iter)> pro_fun = [&](const int & a) {
		return  mvv_pro_mat[from][a] * visited[a];
	};
	for (int idx(0); idx < m_inds.size(); ++idx) {
		auto& cursol = dynamic_cast<TravellingSalesman::SolType&>(*m_inds[idx]);
		std::fill(visited.begin(), visited.end(), 1);
		int first = rnd->uniform.nextNonStd<int>(0, numCity);
		int curIdx(0);
		cursol.variable()[curIdx++] = first;
		visited[first] = 0;
		from = first;
		while (curIdx < numCity) {
			first = rnd->uniform.spinWheel(0, numCity, pro_fun);
			cursol.variable()[curIdx++] = first;
			visited[first] = 0;
			from = first;
		}
	}
	return tag;
}

void ofec::PopGL_NBN_COM_ALG::initialize(Problem* pro, Random* rnd)
{
	Population::initialize(pro, rnd);
	m_rnd = rnd->getSharedPtr();
	mv_dis2parent.resize(m_inds.size());
	std::fill(mv_dis2parent.begin(), mv_dis2parent.end(), 0);
	mv_belong.resize(m_inds.size());
	std::fill(mv_belong.begin(), mv_belong.end(), 0);
	int numCity = CAST_TSP(pro)->numVariables();
	mvv_pro_mat.resize(numCity);
	for (auto& it : mvv_pro_mat) {
		it.resize(numCity);
		std::fill(it.begin(), it.end(), 0);
	}

	mv_sol_edges.resize(m_inds.size());
	m_obj.resize(m_inds.size());
	m_fitness.resize(m_inds.size());
	m_weight.resize(m_inds.size());


}

int ofec::PopGL_NBN_COM_ALG::evaluate(Problem* pro, Algorithm* alg)
{
	auto flag= Population::evaluate(pro, alg);
	for (int idx(0); idx < m_inds.size(); ++idx) {
		CAST_TSP(pro)->transferEdgeSol(*m_inds[idx], mv_sol_edges[idx]);
	}

	udpateProMat(pro);


	std::vector<int> sortSolIds(m_inds.size());
	for (int idx(0); idx < sortSolIds.size(); ++idx) {
		sortSolIds[idx] = idx;
	}

	std::sort(sortSolIds.begin(), sortSolIds.end(), [&](int a,int b) {
		return m_fitness[a] < m_fitness[b];
	});

	std::fill(mv_dis2parent.begin(), mv_dis2parent.end(), std::numeric_limits<double>::max());

	for (int idx(0); idx < sortSolIds.size(); ++idx) {
		int ida = sortSolIds[idx];
		for (int idy(0); idy < idx; ++idy) {
			int idb = sortSolIds[idy];
			double curdis = pro->variableDistance(*m_inds[ida], *m_inds[idb]);
			if (mv_dis2parent[idb] > curdis) {
				mv_dis2parent[idb] = curdis;
				mv_belong[idb] = ida;
			}
			else if (mv_dis2parent[idb] == curdis && m_rnd->uniform.next() < 0.5) {
				mv_belong[idb] = ida;
			}
		}
	}
	
	return flag;
	//return 0;
}

void ofec::PopGL_NBN_COM_ALG::udpateProMat(Problem* pro)
{

	
	{
		double m_memoryMaxObj(0), m_memoryMinObj(0);
		m_memoryMaxObj = m_memoryMinObj = this->m_inds[0]->objective(0);
		for (int i = 0; i < this->size(); ++i) {
			Real obj = this->m_inds[i]->objective(0);
			if (obj > m_memoryMaxObj) m_memoryMaxObj = obj;
			if (obj < m_memoryMinObj) m_memoryMinObj = obj;
		}
		//std::vector<int> indiv(this->m_num_vars);
		Real gap = m_memoryMaxObj - m_memoryMinObj + 1e-5;
		for (int i = 0; i < this->size(); ++i) {
			m_obj[i] = this->m_inds[i]->objective(0);
			if (pro->optMode(0) == OptMode::kMinimize)
				m_fitness[i] = (m_memoryMaxObj - m_obj[i] + 1e-5) / gap;
			else
				m_fitness[i] = (m_obj[i] - m_memoryMinObj + 1e-5) / gap;
			m_weight[i] = 1. / (1 + exp(-m_fitness[i]));
			//m_exMemory[i].push_front(i);
		}



		for (auto& it : mvv_pro_mat) {
			std::fill(it.begin(), it.end(), 0);
		}

		for (int idIndi(0); idIndi < m_inds.size(); ++idIndi) {
			for (int idFrom(0); idFrom < mv_sol_edges[idIndi].size();++idFrom) {
				for (auto& idTo : mv_sol_edges[idIndi][idFrom]) {
					mvv_pro_mat[idFrom][idTo] += m_weight[idIndi];
				}
			}
		}
		double worstWeight = m_weight.front();
		for (auto& it : m_weight) {
			worstWeight = std::min(it, worstWeight);
		}
		for (auto& it : mvv_pro_mat) {
			for (auto& it2 : it) {
				it2 = std::max(it2, worstWeight);
			}
		}

	}
}
