#include "ESs_population.h"
#include <algorithm>

namespace ofec {

	ESs_population::ESs_population(size_t no, size_t dim, size_t size_subpop, Real alpha, Real beta) : population(no), m_probability(CAST_RP_EPANET->number_node()), \
		m_node_data_obj(CAST_RP_EPANET->number_source(), std::make_pair(0.0, 0)), m_size_offspring(size_subpop), m_delta_user(alpha), \
		m_upper_quantile(beta)
	{

	}

	int ESs_population::evolve() {
		int tag = kNormalEval;

		update_probability();
		tag = produce_offspring(m_size_offspring, m_delta_user);
		select(m_upper_quantile);

		++m_iteration;
		return tag;
	}

	int ESs_population::produce_offspring(int size_new, Real detla_user) {
		int tag = kNormalEval;
		m_whole_pop.clear();
		for (size_t i = 0; i < size(); ++i) {
			m_whole_pop.push_back(std::make_pair(m_individuals[i]->variable(), m_individuals[i]->objective()[0]));
		}
		for (size_t i = 0; i < size_new; ++i) {
			int random = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, size());
			m_individuals[random]->mutate(detla_user, m_probability);
			tag = m_individuals[random]->evaluate();
			m_whole_pop.push_back(std::make_pair(m_individuals[random]->variable(), m_individuals[random]->objective()[0]));
			if (tag != kNormalEval) break;
		}
		return tag;
	}

	Real ESs_population::cal_standard_deviation(const std::vector<float>& vec1, const std::vector<float>& vec2) {
		if (vec1.size() != vec2.size()) return -1;
		int size = vec1.size();
		Real sum = 0;
		for (size_t i = 0; i < size; ++i) {
			sum += pow(vec1[i] - vec2[i], 2);
		}
		return sqrt(sum / (Real)size);
	}

	void ESs_population::update_probability() {
		/// update probability of node

		for (auto &i : m_individuals) {
			m_node_data_obj[i->variable().index() - 1].first += i->objective()[0];
			++(m_node_data_obj[i->variable().index() - 1].second);
		}
		std::vector<Real> mean_node(m_node_data_obj.size(), 0.0);
		size_t count_node = 0;
		for (size_t i = 0; i < CAST_RP_EPANET->number_node(); ++i) {
			if (m_node_data_obj[i].second != 0) {
				mean_node[i] = m_node_data_obj[i].first / m_node_data_obj[i].second;
				++count_node;
			}
		}

		/*Real sum = 0.;
		for (auto &i : mean_node)
		sum += i;
		for (size_t i = 0; i < CAST_RP_EPANET->number_node(); ++i) {
		if (count_node > 1 && m_node_data_obj[i].second != 0)
		m_probability[i] = (1 - mean_node[i] / sum) / (count_node - 1);
		else if (count_node == 1 && m_node_data_obj[i].second != 0)
		m_probability[i] = 1;
		else m_probability[i] = 0.001;
		}*/

		Real max = 0;
		for (auto &i : mean_node)
			if (max < i) max = i;
		std::vector<Real> inverse;
		int size_node = mean_node.size();
		for (auto &i : mean_node) {
			if (i == 0 || i == max)
				inverse.push_back(max / size_node);
			else
				inverse.push_back(max - i);
		}
		Real sum = 0;
		for (auto &i : inverse)
			sum += i;
		m_probability.clear();
		for (auto &i : inverse)
			m_probability.push_back(i / sum);
	}

	int ESs_population::select(Real upper_quantile) {
		int tag = kNormalEval;
		std::sort(m_whole_pop.begin(), m_whole_pop.end(), CmpByValue());
		int size_whole = m_whole_pop.size();
		Real critical_point = upper_quantile * size_whole;
		for (size_t i = 0; i < size(); ++i) {
			int random = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, critical_point);
			m_individuals[i]->variable() = m_whole_pop[random].first;
			m_individuals[i]->objective()[0] = m_whole_pop[random].second;
			m_whole_pop[random].first = m_whole_pop[critical_point + i].first;
			m_whole_pop[random].second = m_whole_pop[critical_point + i].second;
		}

		return tag;
	}

}