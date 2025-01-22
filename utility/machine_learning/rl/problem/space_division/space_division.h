#ifndef SPACE_DIVISION_H
#define SPACE_DIVISION_H

#include<cmath>
#include"../problem_base.h"
#include"../../../kd-tree/kdtree_space.h"
#include"../../../../core/global.h"
#include"../../../../core/algorithm/Solution.h"
#include"../../../../core/algorithm/population.h"
#include"../../../../core/problem/continuous/continuous.h"
#include"../../../../utility/parameter/param_map.h"
#include<numeric>
#include<memory>
namespace rl {
	using namespace OFEC;

	enum state2action_mode { global_action_mode, local_action_mode };
	class space_division :public problem_base<int, int> {
		
	public:
		space_division(param_map & v);
		space_division(int state_num = 100);
		virtual ~space_division();
		void initialize()override;
		void set_actions4states()override;//set actions for states;
		void set_state2action_mode(state2action_mode mode) { m_state2action_mode = mode; }
		int get_state(const solution<> & s)const;//return the state attached to the solution;
		const std::vector <std::pair<float, float>> &  space_area(int idx)const;// return the subspace domain;
		float space_volume(int state_idx)const;
		void split_state(int state);//split the state;
		void build_neighborhood();//build the neighborhood (including itself);
		const std::vector<int> & neighbor(int state)const ;//return the neighbor state (including itself);
		int division_depth(int state_idx)const { return m_division_depth[state_idx]; }
		void show()const;
		bool is_in_state(int state,const solution<> & sample)const;
		float direction_deviation(int state_base, int one, int other)const;//calc the cos of two policy direction;
		//void update_state2agent(global * global);
	protected:
		std::vector<float> m_state_ratio; // every status with a division ratio;
		std::vector<std::pair<float, float>>  m_state_space;// all problem search space;
		mutable KDTreeSpace::PartitioningKDTree<float> m_state_division; //to inquire the state of solution;
		std::vector<std::vector<int>> m_neighbor_link;//store the neighbor of every state;
		std::vector<int> m_division_depth;//count the division depth of every state;
		state2action_mode m_state2action_mode;
		//int m_state_num;//the division num of state;
		//rl::td_update<long int, long int> m_rl_update;//the updata alg. of RL;
		//std::map<int, std::vector<int>> m_state2actions;//the state with its possible actions;
		//std::vector<std::vector<std::shared_ptr<OFEC::Solution<>>>> m_state2agent;//the agents in every state;
	private:
		std::vector<int> find_neighbor_state(int state)const;//find the neighbor state;
	};
}

#endif // !SPACE_DIVISION_H

