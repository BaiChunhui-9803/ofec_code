/********* Begin Register Information **********
{
	"name": "SPMOEA2",
	"identifier": "SPMOEA2",
	"problem tags": [ "OOMOP", "ConOP", "MOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 8 June 2023 by Qingshan Tan (email:qingshan.t@cug.edu.cn)
// based on spmoea1_4, just use historical front sols to update space info
// ����ǿ��ѧϰ�ķ�������ӿռ佻��ѡ����Ե����⣬�ӿռ�֮���״̬ת�Ƹ���
// ���ӿռ��ǰ�ؽ�����Ӹ��ʱ���
// ����ȺΪ��λ������ǿ��ѧϰ��������ִ�л���ѡ��
// ���õ���ǿ��ѧϰ
// ��һ��Q�������и����Q��ϲ���ѹ��������״̬�ռ�

#ifndef OFEC_SPMOEA2_H
#define OFEC_SPMOEA2_H

#include "spmoea.h"

namespace ofec {

	class NonStationaryMDP {
	public:
		NonStationaryMDP(int num_states, int num_actions, int horizon)
			: num_states(num_states), num_actions(num_actions), horizon(horizon) {
			generateNonStationaryTransitions();
			generateNonStationaryRewards();
		}

		void generateNonStationaryRewards();
		void generateNonStationaryTransitions();

		std::vector<Real> generateDirichletDistribution(int size);
		Real generateUniformDistribution(Real min, Real max);
		Real generateNormalDistribution(Real mean, Real stddev);

		std::vector<std::vector<std::vector<std::vector<Real>>>> getTransitionProbabilities() const {
			return transition_probabilities;
		}

		std::vector<std::vector<std::vector<Real>>> getActionSelectionProbabilities() const {
			return action_selection_probabilities;
		}

		std::vector<std::vector<std::vector<Real>>> getRewardFunction() const {
			return reward_function;
		}

	private:
		int num_states;
		int num_actions;
		int horizon;
		std::vector<std::vector<std::vector<std::vector<Real>>>> transition_probabilities;
		std::vector<std::vector<std::vector<Real>>> action_selection_probabilities;
		std::vector<std::vector<std::vector<Real>>> reward_function;
	};

	class SPMOEA2 :public SPMOEA {
	public:

		void initialize_();
		void initPop(Problem* pro, Algorithm* alg, Random* rnd);
		void initiObjSpace(Problem* pro);
		void initiVarSpace(Problem* pro);
		void record() override;

		int evolve(Problem* pro, Algorithm* alg, Random* rnd);

		void calFrontBoundRatio(std::vector<Real>& front_bound_ratio);
		void PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro) override;
		int generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, size_t active_type);

		bool spaceSubdivision(Problem* pro, Random* rnd);
		bool spaceSubdivision(size_t space_inx, Problem* pro, Random* rnd);

		std::vector<std::vector<Real>> sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) override;

		void testCoverage(Problem* pro);
		
		void calculateReward(Problem* pro);
		void updateQtable();
		void scaleQtable();
		void updateTransiteProbability();

		void updateInteractiveSols();

		void updateLinkSubspace(size_t inx, const std::vector<size_t>& front_spaces, Random* rnd);

		bool indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol);

		//std::vector<std::pair<size_t,bool>> subspaceSeperable(size_t inx);//����ӿռ��Ƿ���пɷ���
		void splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd);
		//void splitSubspace(size_t inx,size_t num, int flag);
		void splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag);
		int findSplitDim(int inx, Problem* pro);

		void findClusterCenterSsp();

		void updateFrontRegionLinkSpace(Problem* pro, Random* rnd) override;
		std::vector<std::vector<size_t>> linearClusterFrontSpace(std::vector<size_t>& frontspace, Problem* pro, Random* rnd) override;

		void updateSubPopSpace(size_t pop_inx, Problem* pro, Random* rnd);

		//void updateCluster();
		void clusterSubspace();
		std::vector<std::vector<size_t>> clusterFrontSpace(const std::vector<size_t>& frontspace);
		bool subspaceLink(size_t inx1, size_t inx2);

		//void NDSort(Population<Solution<>>& pop);
		void NDSort(std::vector<std::shared_ptr<Solution<>>>& pop);
		//void NDSort(Population<std::unique_ptr<Solution<>>>& pop);
		void recordMetrics(Problem* pro, Algorithm* alg);
		//void recordHisFront(Problem *pro);

		std::vector<Real> getObjUpdateFre() { return m_update_tree; }
		std::vector<std::vector<Real>> getEERatio() { return m_num_e_e; }

		int generateRandomChoice(const std::vector<Real>& probabilities, Random* rnd);

		std::vector<std::vector<Real>> qLearningNonStationaryMDP(const NonStationaryMDP& mdp, Random* rnd,
			Real learning_rate,Real discount_factor,Real exploration_prob,int num_episodes);


#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;

	private:
		std::unique_ptr<SPMOEA_pop> m_rl_pop;
		std::vector<std::unique_ptr<std::map<size_t, std::map<size_t, Real>>>> m_ind_Qtable;//״̬����ֵ����
		//std::unique_ptr<std::map<size_t, std::map<size_t, Real>>> m_immedite_reward;//��ʱ����
		std::vector<std::unique_ptr<std::tuple<size_t, size_t, size_t, Real>>> m_reward;//��ʱ����s-a-s-r��Ԫ��
		bool m_update_cluster = false;
		size_t m_num_series = 10;//���ڻ����仯ʱ������Ԥ��ĳ���
		std::vector<std::pair<Real, Real>> m_pop_range;//��ǰ��Ⱥ��Χ
		std::vector<std::pair<Real, Real>> m_pop_var_range;//��ǰ��Ⱥ��Χ
		std::vector<std::pair<Real, Real>> m_front_pop_range;//ǰ����ȺĿ�귶Χ
		//std::vector <std::shared_ptr<Solution<>>> m_history_front_sols;//store all non-diminated Solutions
		std::vector <std::shared_ptr<Solution<>>> m_archive;//store a number of front sols
		std::vector <std::shared_ptr<Population<Solution<>>>> m_gen_front_sols;//store front sols in each generation

		std::vector<Solution<>> m_subobj_opt_sol;//��Ŀ�����ŵı߽����
		//std::vector<std::unique_ptr<ObjRegionInfo>> m_obj_region_info; //Ŀ���ӿռ���Ϣ
		size_t m_converge_age = 10;//��������������ֵ

		std::vector<size_t> m_over_space_index;//��¼ѡ�������������ӿռ�
		std::vector<std::vector<size_t>> m_front_space_inx;//ǰ���ӿռ�����

		bool m_normalize = true;//Ŀ��ռ��һ�����
		bool m_add_neighbor = false;//�������������ʱ��ǰ���ӿռ�����
		bool m_sample_in_basin = false;//�����ӿռ��ڲ������ǻ������������
		bool m_subspace_select = false;//����ѡ��ʽ
		bool m_evolve_by_potential;//�Ӵ����ɻ���Ǳ���������Ǵ�ͳ��ʽ
		bool m_mean_potential = true;//������Ǳ��Ϊ��ֵ������ֵ
		std::vector<Real> m_update_tree;
		std::vector<std::vector<Real>> m_num_e_e;
		std::vector<size_t> m_his_effect_eval;

		std::vector<std::vector<size_t>> m_region_front_space;
		size_t m_divide_iteration = 0;
		size_t m_divide_granularity;
		size_t m_stage_last_time = 0;

		std::vector<size_t> m_real_front_subspaces;
		std::vector<size_t> m_multi_segment_subspaces;
		std::vector<size_t> m_match_subspace;
		std::vector<size_t> m_error_subspace;
		std::vector<size_t> m_loss_subspace;
		std::vector<size_t> m_match_multi_seg_subspace;

		Real m_learning_rate = 0.1;
		Real m_discount_factor = 0.9;
		Real m_exploration_prob = 0.2;
		int m_num_episodes = 50;

		size_t m_switch_period;
		size_t m_num_rl;
		size_t m_num_extend;

		std::vector<std::vector<size_t>> m_recog_neighs;

		std::vector<std::map<size_t, Real>> m_all_front_nei_prob;//�ӿռ������Ը���


	};

}

#endif  // !OFEC_SPMOEA2_H