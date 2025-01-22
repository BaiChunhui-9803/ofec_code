/******************************************************************************
* Project: Framework for Reinforcement Learning (Also the part of utilities in OFEC)
*******************************************************************************
* Author: Xia Hai
* Email: strawberry9583@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
*  Update_base template serves as all RL Algs. with state_type and aciton type fixed.
*
*********************************************************************************/
#ifndef RL_UPDATE_BASE_H
#define RL_UPDATE_BASE_H

#include"../../../../core/algorithm/algorithm.h"
#include"../model/base.h"
#include"../../parameter/param_map.h"
#include"../../functional.h"
#include<random>
#include<algorithm>
#include<map>
#include<unordered_map>
#include<iostream>
#include<vector>
#include<set>
#include<fstream>

namespace rl {
	struct state_value_type {
	public:
		float getError() { return m_error; }
		float m_state_value;
		float m_error;
		bool operator==(const state_value_type& rhs)const { return m_state_value == rhs.m_state_value&& m_error == rhs.m_error; }
	};

	struct action_value_type {
	public:
		float m_action_value;
		float m_error;
		bool operator==(const action_value_type & rhs)const { return m_action_value == rhs.m_action_value&& m_error == rhs.m_error; }
	};

	struct state_action_pair_hash {
		template<class T1, class T2>
		std::size_t operator() (const std::pair<T1, T2>& p) const {
			auto h1 = std::hash<T1>{}(p.first);
			auto h2 = std::hash<T2>{}(p.second);
			return h1 ^ h2;
		}
	};

	template<typename state_type, typename action_type>
	class update_base :public base {

	public:


	public:
		update_base(const ofec::ParameterMap &v) :base(), m_learning_ratio(v.at("learningRatio")) {
			m_random_engine.seed(0);
		}
		update_base() :base(), m_learning_ratio(0.1) { m_random_engine.seed(0); };
		virtual ~update_base() {};

		//return the states with the biggest v_error_value;
		std::vector<state_type> most_uncertain_state()const {
			std::vector<state_type> state_set;
			float cur_max_error = std::numeric_limits<float>::lowest();
			for (std::pair<state_type, state_value_type>& i : m_v_value) {
				if (cur_max_error < i.second.m_error) {
					state_set.clear();
					state_set.emplace_back(i.first);
					cur_max_error = i.second.m_error;
				}
				else if (cur_max_error == i.second.m_error) {
					state_set.emplace_back(i.first);
				}

			}
			return state_set;
		}

		//return all states with lowest state values from all states;
		std::vector<state_type> worst_state()const {
			std::vector<state_type> all_states(m_v_value.size());
			int idx = 0;
			for (const auto & state : m_v_value) {
				all_states[idx] = state.first;
				++idx;
			}
			return worst_state(all_states);
		}



		//return all states with lowest state values in limited states;
		std::vector<state_type> worst_state(const std::vector<state_type> & states)const {
			std::vector<state_type> worststates;
			float min_statevalue = (std::numeric_limits<float>::max)();
			for (const auto & ele : states) {
				auto value = m_v_value.at(ele).m_state_value;
				if (value < min_statevalue) {
					min_statevalue = value;
					worststates.clear();
					worststates.emplace_back(ele);
				}
				else if (value == min_statevalue) {
					worststates.emplace_back(ele);
				}

			}
			return worststates;

		}

		//randomly choose the worst states with lowest value in limited states;
		state_type random_worst_state(const std::vector<state_type> & states)const {
			std::vector<state_type> all_worst = worst_state(states);
			return random_state(all_worst);
		}

		//return all states with best state values from all states;
		std::vector<state_type> best_state()const {
			std::vector<state_type> all_states(m_v_value.size());
			int idx = 0;
			for (const auto & state : m_v_value) {
				all_states[idx] = state.first;
				++idx;
			}
			return best_state(all_states);
		}


		//all states with best state values in limited states;
		std::vector<state_type> best_state(const std::vector<state_type> & states)const {
			std::vector<state_type> beststates;
			auto max_statevalue = std::numeric_limits<float>::lowest();
			for (const auto & ele : states) {
				auto value = m_v_value.at(ele).m_state_value;
				if (value > max_statevalue) {
					max_statevalue = value;
					beststates.clear();
					beststates.emplace_back(ele);
				}
				else if (m_v_value.at(ele).m_state_value == max_statevalue) {
					beststates.emplace_back(ele);
				}

			}
			return beststates;
		}

		// randomly choose from the limited states;
		state_type random_state(const std::vector<state_type> & states)const {
			std::uniform_int_distribution<int> uni_int(0, (int)states.size() - 1);
			return states[uni_int(m_random_engine)];
		}

		//randomly choose the states with best value in limited states;
		state_type random_best_state(const std::vector<state_type> & states)const {
			std::vector<state_type> all_best = best_state(states);
			return random_state(all_best);
		}

		//In limited states, eplison possibility to choose randomly; 
		//1-epsilon possibility to choose the state with best values;
		state_type epsilon_state(const std::vector<state_type> & states, float epsilon)const {
			std::uniform_real_distribution<float> d_rand;

			if (d_rand(m_random_engine) < epsilon) {
				return random_state(states);
				//uniform_int_distribution<int> uni_int(0, actions.size() - 1);
				//return actions[uni_int(m_random_engine)];
			}
			else {
				return random_best_state(states);
			}

		}

		//greedy action selecting;
		//if there are actions with same values, then return all the actions;
		std::vector<action_type> best_action(const state_type & state) const {
			std::vector<action_type> action_set;
			float cur_max_value = std::numeric_limits<float>::lowest();
			// find all actions with max_action_value;
			for (auto & i : m_q_value) {
				if (i.first.first == state) {
					if (i.second.m_action_value == cur_max_value) {
						action_set.emplace_back(i.second.m_action_value);
					}
					else if (i.second.m_action_value > cur_max_value) {
						action_set.clear();
						action_set.emplace_back(i.first.second);
						cur_max_value = i.second.m_action_value;
					}
				}
			}
			return action_set;
		}

		//greedy action selecting in limited action range;
		std::vector<action_type> best_action(const state_type & state, const std::vector<action_type> & actions)const {
			std::vector<action_type> action_set;
			float cur_max_value = std::numeric_limits<float>::lowest();
			// find all actions with max_action_value;
			for (auto & i : actions) {
				if (cur_max_value == m_q_value.at(std::make_pair(state, i)).m_action_value) {
					action_set.emplace_back(i);
				}
				else if (cur_max_value < m_q_value.at(std::make_pair(state, i)).m_action_value) {
					action_set.clear();
					action_set.emplace_back(i);
					cur_max_value = m_q_value.at(std::make_pair(state, i)).m_action_value;
				}
			}
			if (action_set.empty()) {
				std::cout << "action set empty" << std::endl;
			}
			return action_set;

		}

		action_type random_best_action(const state_type & state) const {

			std::vector<action_type> action_set(best_action(state));
			return random_action(state, action_set);
			//uniform_int_distribution<int> uni_int(0, action_set.size() - 1);
			//return action_set[uni_int(m_random_engine)];
		}

		//greedy action selecting in limited action range;
		action_type random_best_action(const state_type & state, const std::vector<action_type> & actions)const {
			std::vector<action_type> action_set(best_action(state, actions));
			if (action_set.size() == 0) {
				std::cout << "bug occurred" << std::endl;
			}
			return random_action(state, action_set);
			//uniform_int_distribution<int> uni_int(0, action_set.size() - 1);
			//return action_set[uni_int(m_random_engine)];
		}

		//eplison_greedy action selecting; with epsilon probability to choose randomly;
		action_type  epsilon_action(const state_type & state, float epsilon)const {
			//m_epsilon_greedy = epsilon;
			std::vector<action_type> actions;
			for (auto & i : m_q_value) {
				if (i.first.first == state) {
					actions.emplace_back(i.first.second);
				}
			}
			std::uniform_real_distribution<float> d_rand;

			if (d_rand(m_random_engine) < epsilon) {
				std::uniform_int_distribution<int> uni_int(0, (int)actions.size() - 1);
				return actions[uni_int(m_random_engine)];
			}
			else {
				return random_best_action(state);
			}
		}



		//eplison_greedy action selecting; with epsilon probability to choose randomly in limited range;
		action_type  epsilon_action(const state_type & state, float epsilon, const std::vector<action_type> & actions) const {
			//m_epsilon_greedy = epsilon;
			std::uniform_real_distribution<float> d_rand;

			if (d_rand(m_random_engine) < epsilon) {
				return random_action(state, actions);
				//uniform_int_distribution<int> uni_int(0, actions.size() - 1);
				//return actions[uni_int(m_random_engine)];
			}
			else {
				return random_best_action(state, actions);
			}
		}

		//random action selecting;
		action_type random_action(const state_type & state) const {
			std::vector<action_type> actions;
			for (auto & i : m_q_value) {
				if (i.first.first == state) {
					actions.emplace_back(i.first.second);
				}
			}
			std::uniform_int_distribution<int> uni_int(0, (int)actions.size() - 1);
			return actions[uni_int(m_random_engine)];
		}

		//random action selecting in limited range;
		action_type random_action(const state_type & state, const  std::vector<action_type> & actions) const {

			std::uniform_int_distribution<int> uni_int(0, actions.size() - 1);
			//std::map<state_type, std::vector<action_type>>::const_iterator cit = m_state_actions.find(state);
			//if (cit != m_state_actions.end()) {
			//	return cit->second[uni_int(m_random_engine)];
			//	//return m_state_actions.at(state)[uni_int(m_random_engine)];
			//}
			return actions[uni_int(m_random_engine)];
		}


		//select action by probability matching method;
		template<typename T>
		action_type probability_match_action(const state_type & state,
			//			void(*possibility_match)(const std::vector<T> & value, std::vector<double> & possibility)=OFEC::linearly_posibility_match<double> )const {
			std::function<void(const std::vector<T> & value, std::vector<float> & possibility)>possibility_match = OFEC::linearly_possibility_match<float>)const {
			std::vector<action_type> actions;
			std::vector<float> action_values;
			std::vector<float> action_possibility;
			for (auto & i : m_q_value) {
				if (i.first.first == state) {
					actions.emplace_back(i.first.second);
					action_values.emplace_back(i.second.m_action_value);
				}
			}
			/*OFEC::offset_positive(action_values,1.);*/
			possibility_match(action_values, action_possibility);
			//OFEC::linearly_posibility_match(action_values, action_posibility);
			std::uniform_real_distribution<float> ind_dis;
			float cumulative_posi = ind_dis(m_random_engine);
			int action_idx;
			for (action_idx = 0; action_idx < action_values.size(); ++action_idx) {
				cumulative_posi -= action_possibility[action_idx];
				if (cumulative_posi < 0.) {
					break;
				}
			}
			return actions[action_idx];
		}

		//select action by probability matching method in limited range;
		template<typename T>
		action_type probability_match_action(const state_type & state, const  std::vector<action_type> & actions,
			//void(*possibility_match)(const std::vector<T> & value, std::vector<double> & possibility) = OFEC::linearly_posibility_match<double>)const {
			std::function<void(const std::vector<T> & value, std::vector<float> & possibility)>possibility_match = OFEC::linearly_possibility_match<float>)const {
			std::vector<float> action_values;
			std::vector<float> action_possibility;
			for (const auto & _action : actions) {
				action_values.emplace_back(m_q_value.at(std::make_pair(state, _action)).m_action_value);
			}
			//OFEC::offset_positive(action_values,1.);
			//OFEC::normalization(action_values);
			possibility_match(action_values, action_possibility);
			//OFEC::linearly_posibility_match(action_values, action_posibility);
			std::uniform_real_distribution<float> ind_dis;
			float cumulative_posi = ind_dis(m_random_engine);
			int action_idx;
			for (action_idx = 0; action_idx < action_values.size(); ++action_idx) {
				cumulative_posi -= action_possibility[action_idx];
				if (cumulative_posi < 0.) {
					break;
				}
			}
			return actions[action_idx];
		}


		//bool is_state_occur(const state_type & state)const {
		//	//check whether the state_type is;
		//	for (auto & i : m_q_value) {
		//		if (i->first.first == state) {
		//			return true;
		//		}
		//	}
		//	return false;
		//}

		//interface of m_q_value;
		void set_q_value(const std::map<std::pair<state_type, action_type>, action_value_type> & q_value) {
			m_q_value = q_value;
		}
		void set_q_value(const state_type & state, const action_type & action, const action_value_type & q_value) {
			m_q_value[make_pair(state, action)] = q_value;
		}
		const std::unordered_map<std::pair<state_type, action_type>, action_value_type, state_action_pair_hash> & q_value() const { return m_q_value; }
		std::unordered_map<std::pair<state_type, action_type>, action_value_type, state_action_pair_hash>& q_value() { return m_q_value; }
		//interface of m_q_value_error;
		//std::map<std::pair<state_type, action_type>, float>& q_value_error() { return m_q_value_error; }
		//const std::map<std::pair<state_type, action_type>, float> & q_value_error()const { return m_q_value_error; }
		//void set_q_value_error(const std::map<std::pair<state_type, action_type>, float> & q_value_error) {
		//	m_q_value_error = q_value_error;
		//}
		//void set_q_value_error(const state_type & state, const action_type & action, float q_value_error) {
		//	m_q_value_error[make_pair(state, action)] = q_value_error;
		//}
		//interface of m_v_value;
		std::unordered_map<state_type, state_value_type> & v_value() { return m_v_value; }
		const std::unordered_map<state_type, state_value_type>& v_value()const { return m_v_value; }
		void set_v_value(const std::map<state_type, state_value_type> & v_value) { m_v_value = v_value; }
		void set_v_value(const state_type & state, const state_value_type & v_value) {
			m_v_value[state] = v_value;
		}
		//interface of m_v_value_error;
		//std::map<state_type, float> & v_value_error() { return m_v_value_error; }
		//const std::map<state_type, float> & v_value_error()const { return m_v_value_error; }
		//void set_v_value_error(const std::map<state_type, float> & v_value_error) { m_v_value_error = v_value_error };
		//void set_v_value_error(const state_type & state, float v_value_error) {
		//	m_v_value_error[state] = v_value_error;
		//}
		//interface of m_state_actions;
		//void set_state_actions(const map <state_type, vector<action_type> > & state2actions) {
		//	m_state_actions = state2actions;
		//}

		//interface of m_state2actions;
		//const map<state_type, vector<action_type>> & state2actions() const { return m_state2actions; }
		//map<state_type, vector<action_type>>& state2actions() { return m_state2actions; }

		//interface of m_learning_ratio;
		void set_learning_ratio(float ratio) { m_learning_ratio = ratio; }
		float learning_ratio()  const { return m_learning_ratio; }

		//interface of m_epsilon_greedy;
		//void set_epsilon(float epsilon) { m_epsilon_greedy = epsilon; }
		//float  epsilon_greedy() { return m_epsilon_greedy; }

		//initialize the state value, action value and m_state_actions;
		//@param: state_actions: actions can be taken in specific states; 
		//		  initial_v: state values; initial_v_error: state value errors;
		//        initial_q: state action values; initila_q_error: state action values errors;
		//		  if_set_actions: whether set action values for states; (save memory when the number of states is large)
		//						  flase--> not set action values; true--> set action values; 
		void initialize(const std::unordered_map<state_type, std::vector<action_type> >& state_actions, \
			float initial_v = 0.f, float initial_v_error = 0.f, float initial_q = 0.f, float initial_q_error = 0.f, \
			bool if_set_actions = true) {
			//std::cout << "update_base initialize begin" << std::endl;
			for (auto i = state_actions.begin(); i != state_actions.end(); ++i) {
				m_v_value[i->first].m_state_value = initial_v;
				m_v_value[i->first].m_error = initial_v_error;
				if (if_set_actions) {
					//the initial estimation error for v_value is infinite;
					//m_v_value_error[i->first] = std::numeric_limits<float>::max();
					for (auto & j : i->second) {
						m_q_value[std::make_pair(i->first, j)].m_action_value = initial_q;
						m_q_value[std::make_pair(i->first, j)].m_error = initial_q_error;

						//the initial estimation error for q_value is infinite;
						//m_q_value_error[std::make_pair(i->first, j)] = std::numeric_limits<float>::max();
					}

				}
			}

			//std::cout << "update_base initialize ends" <<m_i++<< std::endl;
		}

		//float action_entropy() {
		//	for(auto & i:m_)
		//}
	protected:


	protected:
		//the state and its conresponding actions;
		//map<state_type, vector<action_type> > m_state_actions;
		//the state_value;
		std::unordered_map<state_type, state_value_type> m_v_value;
		//the state_action value;
		std::unordered_map<std::pair<state_type, action_type>, action_value_type, state_action_pair_hash> m_q_value;
		//the action_value error;
		//std::map<std::pair<state_type, action_type>, float> m_q_value_error;

		//learning ratio;
		float m_learning_ratio;
		//mutable (for the const function access safety) random engine for learning alg.;
		mutable std::default_random_engine m_random_engine;
		static int m_i;
		//epsilon greedy probability;
		//float m_epsilon_greedy;
	};

	template<typename state_type, typename action_type>
	int update_base<state_type, action_type>::m_i = 1;
}
#endif // !UPDATE_ALG_H
