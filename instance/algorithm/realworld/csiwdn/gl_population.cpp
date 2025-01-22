#include "gl_population.h"
#include "../../../problem/realworld/csiwdn/csiwdn.h"

namespace ofec {

	GLPopulation::GLPopulation(size_t no, Environment* env, size_t dim) :
		Population(no, env, dim) {
		auto pro = env->problem();
		m_probability.resize(CAST_CSIWDN(pro)->numberSource(), std::vector<Real>(CAST_CSIWDN(pro)->numberNode()));
		//m_offspring(no, std::vector<int> (CAST_CSIWDN(pro)->numberSource())),
		m_node_data_obj.resize(CAST_CSIWDN(pro)->numberSource(), std::vector<std::pair<Real, size_t>>(CAST_CSIWDN(pro)->numberNode(), std::make_pair(0.0, 0)));
	}

	void GLPopulation::initialize(Environment *env, Random *rnd) {
		Population::initialize(env, rnd);
		auto pro = env->problem();

		m_node_cluster = CAST_CSIWDN(pro)->currentCluster();
	}

	void GLPopulation::evolve(Environment *env, Random *rnd,  bool is_stable, const std::pair<int, int>& source_index) {
		 updateProbability(env, source_index);
		 mutate(env, rnd, source_index);
		 select(env, false, source_index);
	 }

	 void  GLPopulation::mutate(Environment *env, Random *rnd, const std::pair<int, int>& source_index) {
		 auto pro = env->problem();
		 int q = source_index.first, z = source_index.second;
		 for (int j = q; j <= z; ++j) {
			 for (int i = 0; i < m_individuals.size(); i++) {
				 int start = m_individuals[i]->variable().index(j) - 1;
				 int end = *rnd->uniform.nextElem(m_node_cluster.begin(), m_node_cluster.end());
				 if (start == end)
					 m_individuals[i]->trial().variable().index(j) = start + 1;
				 else {
					 auto path_node = CAST_CSIWDN(pro)->shortestPath(start, end);
					 std::vector<Real> p(path_node.size());
					 Real delt = (m_probability[j][end] - m_probability[j][start]) / (path_node.size() - 1);
					 for (size_t k = 0; k < p.size(); k++) {
						 p[k] = m_probability[j][start] + delt * k;
					 }
					 m_individuals[i]->trial().variable().index(j) = path_node[rnd->uniform.spinWheel(p.begin(), p.end()) - p.begin()] + 1;
				 }
			 }
		 }
		 //for (int j = q; j <= z; ++j) {
			// for (int i = 0; i < m_offspring.size(); i++) {
			//	 m_individuals[i]->variable().index(j) = m_offspring[i][j] + 1;
			// }
		 //}

		//for (size_t i = 0; i < size(); ++i) {
		//	m_individuals[i]->mutateFirstPart(m_probability,pro,rnd, source_index);
		//}
	}

	void  GLPopulation::updateProbability(Environment *env, const std::pair<int, int>& source_index) {

		auto pro = env->problem();
		/// update probability of node
		int z = source_index.second, q = source_index.first;
		for (auto &i : m_individuals) {
			for (int j = q; j <= z; j++) {
				if (CAST_CSIWDN(pro)->indexFlag(i->variable().index(z) - 1)) {
					m_node_data_obj[j][i->variable().index(j) - 1].first += i->objective()[0];
					++(m_node_data_obj[j][i->variable().index(j) - 1].second);
				}
			}
		
		}

		std::vector<std::vector<Real>> mean_node(CAST_CSIWDN(pro)->numberSource(), std::vector<Real>(CAST_CSIWDN(pro)->numberNode(), 0.0));
		for (size_t i = q; i <= z; ++i) {
			for (size_t j = 0; j < CAST_CSIWDN(pro)->numberNode(); ++j) {
				if (m_node_data_obj[i][j].second != 0)
					mean_node[i][j] = m_node_data_obj[i][j].first / m_node_data_obj[i][j].second;
				else
					mean_node[i][j] = -1;
			}
		}

		m_probability.clear();
		m_probability.resize(mean_node.size());
		Real max, min;
		for (int j = q; j <= z; j++) {
			max = 0, min = mean_node[j][0];
			for (auto& i : mean_node[j]) {
				if (max < i) max = i;
				if (min > i) min = i;
			}
			std::vector<Real> inverse;
			int size_node = mean_node[j].size();
			for (auto& i : mean_node[j]) {
				if (i == -1)
					i = max;
				inverse.push_back((max - i) / (max - min) + 0.01);
			}

			Real sum = 0;
			for (auto& i : inverse)
				sum += i;

			for (auto& i : inverse)
				m_probability[j].push_back(i / sum);
		}
	}

	void  GLPopulation::fillSolution(VarCSIWDN &indi, Environment *env, const std::pair<int, int>& source_index) {
		for (auto &i : m_individuals)
			i->coverSecondPart(indi, env, source_index);
	}

	void  GLPopulation::select(Environment *env, bool is_stable, const std::pair<int, int>& source_index) {
		for (auto &i : m_individuals)
			i->select(env,is_stable, source_index);

	}

	bool  GLPopulation::isFeasiblePopulation(Environment *env,const Real tar) {
		auto pro = env->problem();
		size_t phase =  CAST_CSIWDN(pro)->phase();
		size_t interval =  CAST_CSIWDN(pro)->interval();
		size_t num = phase*interval;
		Real temp = 0;

		for (size_t i = 0; i < num; ++i) {
			for (int j = 0; j <  CAST_CSIWDN(pro)->numSensor(); ++j) {
				temp += pow( CAST_CSIWDN(pro)->observationConcentration()[j][i], 2);
			}
		}

		Real benchmark = tar * sqrt(temp / ( CAST_CSIWDN(pro)->numSensor()*num));

		size_t count_feasible = 0, count_infeasible = 0;
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			if (m_individuals[i]->objective()[0] <= benchmark) ++count_feasible;
			else ++count_infeasible;
		}

		return count_feasible >= count_infeasible;


	}
}