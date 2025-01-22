#include "gl_calculator.h"


void ofec::GL_calculator::initialize(ofec::Problem* pro) {
	mvv_pro_mat.resize(pro->numberVariables());
	for (auto& it : mvv_pro_mat) {
		it.resize(pro->numberVariables());
	}
}

void ofec::GL_calculator::calWeight(
	std::vector<double>& weight, 
	const std::vector<double>& objs, 
	ofec::Problem* pro) {

	weight.resize(objs.size());
	double m_memoryMaxObj(0), m_memoryMinObj(0);
	m_memoryMaxObj = m_memoryMinObj = objs.front();
	for (int i = 0; i < objs.size(); ++i) {
		Real obj = objs[i];
		if (obj > m_memoryMaxObj) m_memoryMaxObj = obj;
		if (obj < m_memoryMinObj) m_memoryMinObj = obj;
	}
	//std::vector<int> indiv(this->m_num_vars);
	double gap = m_memoryMaxObj - m_memoryMinObj + 1e-5;
	double fit(0);
	for (int i = 0; i < objs.size(); ++i) {
		if (pro->optimizeMode(0) == OptimizeMode::kMinimize)
			fit= (m_memoryMaxObj - objs[i] + 1e-5) / gap;
		else
		fit= (objs[i] - m_memoryMinObj + 1e-5) / gap;
		weight[i] = 1. / (1 + exp(-fit));
		//m_exMemory[i].push_front(i);
	}


}



void ofec::GL_calculator::udpateProMat(
	const std::vector<std::vector<std::vector<int>>*>& edges,
	std::vector<double>& weight) {
	

	for (auto& it : mvv_pro_mat) {
		std::fill(it.begin(), it.end(), 0);
	}

	for (int idIndi(0); idIndi < weight.size(); ++idIndi) {
		auto& curedge = *edges[idIndi];
		for (int idFrom(0); idFrom < curedge[idIndi].size(); ++idFrom) {
			for (auto& idTo : curedge[idFrom]) {
				mvv_pro_mat[idFrom][idTo] += weight[idIndi];
			}
		}
	}
	
}


void ofec::GL_calculator::generateTSPsolution(
	std::vector<int>& tour, 
	ofec::Problem* pro,
	ofec::Random* rnd
	) {
	int numCity = pro->numberVariables();
	std::vector<int> visited(numCity);
	//	std::vector<double> pro;
	int from = 0;
	std::function< double(const int& cur_iter)> pro_fun = [&](const int& a) {
		return  mvv_pro_mat[from][a] * visited[a];
	};
	std::function< double(const int& cur_iter)> pro_fun2 = [&](const int& a) {
	   return visited[a];
	};

	tour.resize(numCity);
	std::fill(visited.begin(), visited.end(), 1);
	int first = rnd->uniform.nextNonStd<int>(0, numCity);
	int curIdx(0);
	tour[curIdx++] = first;
	visited[first] = 0;
	from = first;
	while (curIdx < numCity) {
		first = rnd->uniform.spinWheel(0, numCity, pro_fun);
		if (first == numCity) {
			first = rnd->uniform.spinWheel(0, numCity, pro_fun2);
		}
		tour[curIdx++] = first;
		visited[first] = 0;
		from = first;
	}

}

void ofec::GL_calculator::generateTSPsolution(
	std::vector<int>& tour, 
	ofec::Problem* pro, ofec::Random* rnd, 
	const std::vector<std::vector<int>>& centerSol, double maxRadius)
{
    

}
