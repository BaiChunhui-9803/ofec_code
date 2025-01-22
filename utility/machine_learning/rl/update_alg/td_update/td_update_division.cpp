//#include"td_update_division.h"
//
//namespace rl {
//	td_update_division::td_update_division(OFEC::param_map & v): td_update<int,int>(v) {
//	}
//
//	td_update_division::~td_update_division() {
//	}
//
//	const std::unordered_map<std::pair<int, int>, int,rl::state_action_pair_hash>& td_update_division::visit_count() const {
//		// TODO: 在此处插入 return 语句
//		return m_visit_count;
//	}
//
//	std::unordered_map < std::pair<int, int>, int, rl::state_action_pair_hash > & td_update_division::visit_count() {
//		// TODO: 在此处插入 return 语句
//		return m_visit_count;
//	}
//
//	void td_update_division::initialize(const std::unordered_map<int, std::vector<int>>& state_actions,
//		float initial_v, float initial_v_error, float initial_q, float initial_q_error ){
//		//initialize the visit_count;
//		//for (const auto & i : state_actions) {
//		//	for (const auto & j : i.second) {
//		//		m_visit_count[make_pair(i.first,j)] = 0;
//		//	}
//		//}
//		td_update::initialize(state_actions,initial_v,initial_v_error,initial_q,initial_q_error);
//	}
//
//
//	void td_update_division::add_action_value(std::unordered_map<std::pair<int, int>, rl::action_value_type,rl::state_action_pair_hash> & added_value) {
//		for (auto & i : added_value) {
//			m_q_value[i.first] = i.second;
//		}
//	}
//
//	
//	void td_update_division::splite_state(int state,const std::vector<int> & neighbor_link) {
//		//add the state value;
//		m_v_value[m_v_value.size()] = m_v_value[state];
//		//every state to the new added state;
//		int state_num = m_v_value.size();
//		for (int i = 0; i < state_num; ++i) {
//			if (i == state_num-1) {
//				m_q_value[std::make_pair(state_num - 1, state_num - 1)] = m_q_value[std::make_pair(state, state)];
//			}
//			else {
//				m_q_value[std::make_pair(i,state_num-1)] = m_q_value[std::make_pair(i, state)];
//			}
//		}
//		//new added state to every neighbor states;
//		for (const auto & i : neighbor_link) {
//			m_q_value[std::make_pair(state_num - 1, i)] = m_q_value[std::make_pair(state, i)];
//		}
//	}
//	void td_update_division::update_q(const int & state, const int & action, const int & next_state, float r) {
//		td_update::update_q(state, action, next_state, r);
//		++m_visit_count[std::make_pair(state, action)];
//	}
//}