#include "state_transition.h"
#include "../linear_algebra/eigen.h"

namespace ofec::population_distribution {
	StateTransition::StateTransition(size_t number_states) :
		m_number_states(number_states) {}

	size_t StateTransition::numberStates() const {
		return m_number_states;
	}

	void StateTransition::record(size_t id_from, size_t id_to) {
		if (m_count[id_from].count(id_to)) {
			m_count[id_from][id_to]++;
		}
		else {
			m_count[id_from][id_to] = 1;
		}
		if (!m_sampled_states_ids.count(id_from)) {
			m_state_id_to_order.emplace(id_from, m_order_to_state_id.size());
			m_order_to_state_id.emplace_back(id_from);
			m_sampled_states_ids.insert(id_from);
		}
		if (!m_sampled_states_ids.count(id_to)) {
			m_state_id_to_order.emplace(id_to, m_order_to_state_id.size());
			m_order_to_state_id.emplace_back(id_to);
			m_sampled_states_ids.insert(id_to);
		}
	}

	size_t StateTransition::count(size_t id_from, size_t id_to) const {
		if (m_count.count(id_from) && m_count.at(id_from).count(id_to)) {
			return m_count.at(id_from).at(id_to);
		}
		else {
			return 0;
		}
	}

	void StateTransition::calculate() {
		m_probability.resize(m_state_id_to_order.size());
		for (size_t i = 0; i < m_state_id_to_order.size(); i++) {
			m_probability[i].resize(m_state_id_to_order.size());
			Real sum_count = 0;
			for (size_t j = 0; j < m_state_id_to_order.size(); j++) {
				sum_count += count(m_order_to_state_id[i], m_order_to_state_id[j]);
			}
			if (sum_count > 0) {
				for (size_t j = 0; j < m_state_id_to_order.size(); j++) {
					m_probability[i][j] = count(m_order_to_state_id[i], m_order_to_state_id[j]) / sum_count;
				}
			}
			else {
				for (size_t j = 0; j < m_state_id_to_order.size(); j++) {
					m_probability[i][j] = 1.0 / m_number_states;
				}
			}
		}
	}

	Real StateTransition::probability(size_t id_from, size_t id_to) const {
		if (m_sampled_states_ids.count(id_from)) {
			if (m_sampled_states_ids.count(id_to)) {
				return m_probability[m_state_id_to_order.at(id_from)][m_state_id_to_order.at(id_to)];
			}
			else {
				return 0.0;
			}
		}
		else {
			return 1.0 / m_number_states;
		}
	}

	const std::vector<size_t>& StateTransition::sampledStatesIDs() const {
		 return m_order_to_state_id;
	}

	size_t StateTransition::stateIDtoOrder(size_t state_id) const {
		return m_state_id_to_order.at(state_id);
	}

	std::map<size_t, Real> StateTransition::convergenceProbability(size_t state_id) const {
		Real epsilon = 1e-4;
		Eigen::RowVectorXd current_probability(m_number_states);
		current_probability.setZero();
		current_probability(state_id) = 1;
		Eigen::MatrixXd transition_probability(m_number_states, m_number_states);
		for (size_t i = 0; i < m_number_states; ++i) {
			for (size_t j = 0; j < m_number_states; ++j) {
				transition_probability(i, j) = probability(i, j);
			}
		}
		while(true) {
			auto new_probability = current_probability * transition_probability;
			Real max_error = (new_probability - current_probability).cwiseAbs().maxCoeff();
			//std::cout << max_error << std::endl;
			if (max_error < epsilon) {
				break;
			}
			current_probability = new_probability;
		}
		std::map<size_t, Real> result;
		for (size_t i = 0; i < m_number_states; ++i) {
			result[i] = current_probability(i);
		}
		return result;
	}
}
