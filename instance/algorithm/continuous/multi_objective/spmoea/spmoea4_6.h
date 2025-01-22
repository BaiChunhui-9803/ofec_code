/********* Begin Register Information **********
{
	"name": "SPMOEA4_6",
	"identifier": "SPMOEA4_6",
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
// ��spmoea4_4�Ļ����ϣ��ҵ����ѵ�����󣬷�����Դ�����������ݻ�
// ����Ⱥ�ݻ�ʱ�����оֲ�����

#ifndef OFEC_SPMOEA4_6_H
#define OFEC_SPMOEA4_6_H

#include "spmoea.h"
#include<queue>

namespace ofec {

	class SPMOEA4_6 :public SPMOEA {
	public:

		void initialize_();
		void initPop(Problem* pro, Algorithm* alg, Random* rnd);
		void initiObjSpace(Problem* pro);
		void initiVarSpace(Problem* pro);
		void record() override;

		int evolve(Problem* pro, Algorithm* alg, Random* rnd);
		void assignPop(Population<Solution<>>& pop, Problem* pro, Algorithm* alg, Random* rnd);

		size_t selectSubspace(Problem* pro, Random* rnd);
		std::map<Real, std::pair<size_t, size_t>> selectSubspaceFromObj(Problem* pro, Random* rnd);

		void calFrontBoundRatio(std::vector<Real>& front_bound_ratio);
		void PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro);

		void generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, const std::vector<size_t>& pop_resource, std::vector<int> type);
		std::vector<std::vector<Real>> sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleLinearInFrontSpace(size_t space, size_t ind_inx, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> samplePush(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleByIndPos(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleByPush(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd);

		bool spaceSubdivision(Problem* pro, Random* rnd);

		std::vector<std::vector<Real>> sampleDiverse(Population<Solution<>>& pop, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);

		std::vector<std::vector<Real>> sampleRandom(Population<Solution<>>& pop, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);

		bool indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol);

		void adaptiveSplit(Problem* pro);
		void adaptiveSplitFrontSpace(Problem* pro, Random* rnd);
		//std::vector<std::pair<size_t,bool>> subspaceSeperable(size_t inx);//����ӿռ��Ƿ���пɷ���
		void splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd);
		//void splitSubspace(size_t inx,size_t num, int flag);
		void splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag);
		int findSplitDim(int inx, Problem* pro);

		void findClusterCenterSsp();
		void updateLinkSubspace(size_t inx, const std::vector<size_t>& front_spaces, Random* rnd);

		void testCoverage(Problem* pro);
		void updateNeighSpace();

		void SolutionSelection(Problem* pro, Random* rnd) override;

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


#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;

	private:
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
		size_t m_switch_period;
		size_t m_num_push;
		size_t m_num_extend;
		size_t m_neighs = 5;
		size_t m_immigrant_times = 0;

		std::vector<size_t> m_real_front_subspaces;
		std::vector<size_t> m_real_all_front_subspaces;
		std::vector<size_t> m_multi_segment_subspaces;
		std::vector<size_t> m_match_subspace;
		std::vector<size_t> m_error_subspace;
		std::vector<size_t> m_loss_subspace;
		std::vector<size_t> m_match_multi_seg_subspace;

		std::vector<std::vector<size_t>> m_recog_neighs;

		std::vector<std::map<size_t, Real>> m_all_front_nei_prob;//�ӿռ������Ը���

		size_t m_total_push_times = 0;
		std::vector<size_t> m_push_results = { 0,0,0 };//dominant,dominated,nondominated

		size_t m_total_extend_times = 0;
		std::vector<size_t> m_extend_results = { 0,0,0 };

		std::vector<size_t> m_total_results = { 0,0,0 };//�ۻ��ݻ�ͳ��

		size_t m_push_behind_empty = 0;

		size_t m_queue_size = 20;

		Real m_explore_rate;

		size_t m_explore_fail_times = 0;

		bool m_explore_reset_flag = false;

		std::vector<std::pair<Real, Real>> m_reset_box;

		std::map<size_t, std::vector<size_t>> m_neigh_spaces;

		std::vector<size_t> m_ind_operator_times = { 0,0,0,0,0,0,0,0 };
		std::vector<std::vector<size_t>> m_operator_results;

		std::map<size_t, std::vector<size_t>> m_subspace_pop_index;

		std::vector<std::vector<size_t>> m_pop_index_clusters;
		std::vector<std::vector<size_t>> m_pop_index;

		std::map<size_t, std::pair<size_t, std::queue<int>>> m_subspace_succ_rate;//space_idx,count,succ_rate

		std::map<size_t, std::pair<size_t, std::queue<int>>> m_subspace_interactive_result;//space_idx,count,succ_rate
	};

}

#endif  // !OFEC_SPMOEA4_6_H