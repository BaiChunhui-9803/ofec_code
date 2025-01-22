#ifndef OFEC_POPULATION_DISTRIBUTION_STATE_TRANSITION_H
#define OFEC_POPULATION_DISTRIBUTION_STATE_TRANSITION_H

#include <vector>
#include <map>
#include "../../core/definition.h"
#include <set>

namespace ofec::population_distribution {
	class StateTransition {
	private:
		std::map<size_t, std::map<size_t, size_t>> m_count;		// [i][j] denotes the times state change from i to j
		std::vector<std::vector<Real>> m_probability;
		size_t m_number_states;
		std::vector<size_t> m_order_to_state_id;
		std::set<size_t> m_sampled_states_ids;
		std::map<size_t, size_t> m_state_id_to_order;
	public:
		StateTransition(size_t number_states);
		size_t numberStates() const;
		void record(size_t id_from, size_t id_to);
		size_t count(size_t id_from, size_t id_to) const;
		void calculate();
		Real probability(size_t id_from, size_t id_to) const;
		const std::vector<size_t>& sampledStatesIDs() const;
		size_t stateIDtoOrder(size_t state_id) const;
		std::map<size_t, Real> convergenceProbability(size_t state_id) const;
	};
}

#endif // !OFEC_POPULATION_DISTRIBUTION_STATE_TRANSITION_H
