#include "space_division.h"

namespace rl {
	space_division::space_division(param_map & v):problem_base<int, int>(v.at("stateNum")), \
		m_state_division(CONTINUOUS_CAST->variable_size(), m_state_ratio, m_state_space),\
		m_division_depth(v.at("stateNum"),0) {
	}
	space_division::space_division(int state_num):problem_base<int,int>(state_num),\
		m_state_division(CONTINUOUS_CAST->variable_size(), m_state_ratio, m_state_space) {
	}

	space_division::~space_division() {
	}

	void space_division::initialize() {
		//to allocate ratio for every state evenly;
		for (int i = 0; i < m_state_num; ++i) {
			m_state_ratio.emplace_back(1. / m_state_num);
		}
		// to get the search spaces;
		for (size_t i = 0; i < global::ms_global->m_problem->variable_size(); i++) {
			m_state_space.emplace_back(CONTINUOUS_CAST->range(i).first,
				CONTINUOUS_CAST->range(i).second);
		}
		m_state_division.inputData(m_state_ratio);
		m_state_division.setInitBox(m_state_space);
		//build the k-d state tree, then can be equired.
		m_state_division.buildIndex();
		//m_state_division.regionShow();

		//set states set and action set;
		m_states.resize(m_state_num);
		m_actions.resize(m_state_num);
		std::iota(m_states.begin(), m_states.end(), 0);
		std::iota(m_actions.begin(), m_actions.end(), 0);

		//build the neighborhoood;
		build_neighborhood();

		//set states to actions(set actions for states)
		set_actions4states();

		float neighobr_num = 0;
		for (auto & neigho : m_neighbor_link) {
			neighobr_num += neigho.size();
		}
		
		std::cout << "average No. of neighbors: " << neighobr_num / m_state_num << std::endl;
		std::cout << "neighbor ratio: " << neighobr_num / (m_state_num*m_state_num) << std::endl;
		//std::cin.get();
	}

	void space_division::set_actions4states() {
		//set to ations for states;
		if (m_state2action_mode == state2action_mode::global_action_mode) {
			for (int i = 0; i < m_state_num; ++i) {
				std::vector<int> temp_actions(m_state_num);
				std::iota(temp_actions.begin(), temp_actions.end(), 0);
				this->m_states2actions.emplace(i, temp_actions);
			}
		}
		else if(m_state2action_mode==state2action_mode::local_action_mode) {
			for (int i = 0; i < m_state_num; ++i) {
				this->m_states2actions.emplace(i, m_neighbor_link[i]);
			}
		}
	}

	int space_division::get_state(const solution<>& s) const{
		return m_state_division.get_regionIdx(s.variable().vect());
	}
	const std::vector <std::pair<float, float>> & space_division::space_area(int idx) const {
		return m_state_division.get_box(idx);
	}
	float space_division::space_volume(int state_idx) const {
		return m_state_division.getBoxVolume(state_idx);
	}
	std::vector<int> space_division::find_neighbor_state(int state) const {
		std::vector<int> neighbor;
		m_state_division.find_neighbor(state, neighbor);
		return neighbor;
	}

	void space_division::split_state(int state) {
		//split the volume space in two half;
		if (state >= 0 && state < m_state_num) {
			m_state_division.split_region(state);
			//new state and action added;
			++m_state_num;
			m_states.emplace_back(m_state_num-1);
			m_states2actions.emplace(m_state_num-1, m_actions);
			m_actions.emplace_back(m_state_num-1);
			//update the state2actions;
			for (auto & i : m_states2actions) {
				i.second.emplace_back(m_state_num-1);
			}
			//update division ratio;
			m_state_ratio[state] /= 2.;
			m_state_ratio.emplace_back(m_state_ratio[state]);
		}
		else {
			std::cout << "the splited state not existing" << std::endl;
		}
		//update the division depth of the two halves;
		++m_division_depth[state];
		m_division_depth.emplace_back(m_division_depth[state]);

		//update the neighborhood;
		for (auto & i : m_neighbor_link[state]) {
			if (i != state) {
				std::vector<int> neighbor;
				m_state_division.find_neighbor(i, neighbor);
				neighbor.emplace_back(i);
				m_neighbor_link[i] = neighbor;
			}
		}
		//add the neighborhood of the first half state;
		std::vector<int> neighbor_1;
		m_state_division.find_neighbor(state, neighbor_1);
		neighbor_1.emplace_back(state);
		m_neighbor_link[state] = neighbor_1;
		//add the neighborhood of the second half state;
		std::vector<int> neighbor_2;
		m_state_division.find_neighbor(m_state_num-1, neighbor_2);
		neighbor_2.emplace_back(m_state_num-1);
		m_neighbor_link.emplace_back(neighbor_2);
	}

	void space_division::build_neighborhood() {
		for (int i = 0; i < m_state_num; ++i) {
			auto neighbor_i = find_neighbor_state(i);
			//push the self idx link to the neighbor;
			neighbor_i.push_back(i);
			m_neighbor_link.emplace_back(neighbor_i);
		}
	}
	const std::vector<int>& space_division::neighbor(int state)const {
		return m_neighbor_link[state];
	}

	void space_division::show() const {
		m_state_division.regionShow();
	}

	bool space_division::is_in_state(int state,const solution<> & sample) const {
		const auto & domain =m_state_division.get_box(state);
		for (size_t dim = 0; dim < domain.size(); ++dim) {
			if (sample.variable()[dim] < domain[dim].first || sample.variable()[dim] > domain[dim].second) {
				return false;
			}
		}
		return true;
	}

	float space_division::direction_deviation(int state_base, int one, int other) const {
		const auto base_area = space_area(state_base);
		const auto one_area = space_area(one);
		const auto other_area = space_area(other);
		std::pair<float, float> v_base{ (base_area[0].first + base_area[0].second)*0.5f,
			(base_area[1].first + base_area[1].second)*0.5f };
		std::pair<float,float> v_1{ (one_area[0].first + one_area[0].second)*0.5f,
			(one_area[1].first + one_area[1].second)*0.5f };
		std::pair<float,float>  v_2{ (other_area[0].first + other_area[0].second)*0.5f,
			(other_area[1].first + other_area[1].second)*0.5f };
		v_1.first -= v_base.first;
		v_1.second -= v_base.second;
		v_2.first -= v_base.first;
		v_2.second -= v_base.second;
		float mode_1 = std::sqrt(v_1.first*v_1.first + v_1.second*v_1.second);
		float mode_2 = std::sqrt(v_2.first*v_2.first + v_2.second*v_2.second);
		if (mode_1 == 0.f&&mode_2 == 0.f) {
			return 1.f;
		}
		else if (mode_1 == 0.f || mode_2 == 0.f) {
			return std::cos((2 * OFEC_PI) / neighbor(state_base).size());
		}
		else {
			return (v_1.first*v_2.first + v_1.second*v_2.second) / (mode_1*mode_2);
		}
	}

	//void space_division::update_state2agent(global * global) {
	//	//erase the old mapping from state to agents;
	//	for (auto & i : m_state2agent)
	//		i.clear();
	//	//building the new mapping;
	//	for (auto & i : dynamic_cast<*OFEC::population<>>())
	//		m_state2agent[i->cur_state()].push_back(std::make_shared<Solution<>>(*i));
	//}
}
