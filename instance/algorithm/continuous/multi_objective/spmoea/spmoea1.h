/********* Begin Register Information **********
{
	"name": "SPMOEA1",
	"identifier": "SPMOEA1",
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
// 在spmoea4_8的基础上，根据目标空间和搜索空间解的分布，不平衡分配搜索资源
// 取消额外的精搜子种群，设计贪心概率，在识别的区域精搜
// 根据前沿解的数量或者IGD的变化调节环境选择和贪心概率
// 根据解集的指标变化来触发精搜

// 优化自适应分割子空间的方式

#ifndef OFEC_SPMOEA1_H
#define OFEC_SPMOEA1_H

#include "spmoea.h"
#include "spaceGraph.h"
#include<queue>

namespace ofec {

	class SPMOEA1 :public SPMOEA {
	public:

		void initialize_();
		void initPop(Problem* pro, Algorithm* alg, Random* rnd);
		void initiObjSpace(Problem* pro);
		void initiVarSpace(Problem* pro);
		void record() override;

		int evolve(Problem* pro, Algorithm* alg, Random* rnd);
		void assignPop(Population<Solution<>>& pop, Problem* pro, Algorithm* alg, Random* rnd);

		size_t selectSubspace(Problem* pro, Random* rnd);
		std::map<Real, std::pair<size_t, size_t>> selectSubspaceFromObj(Problem* pro, Random* rnd, int source);
		std::vector<size_t> selectSubspaceFromVar(Problem* pro, Random* rnd);

		void clusterInVarspace(Problem* pro);

		void calFrontBoundRatio(std::vector<Real>& front_bound_ratio);
		void PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro);

		void uniformSample(size_t pop_inx, Problem* pro, Algorithm* alg, Random* rnd);
		void biasSample(size_t pop_inx, Problem* pro, Algorithm* alg, Random* rnd, Real prob);

		void generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, const std::vector<size_t>& pop_resource);
		std::vector<std::vector<Real>> sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleLinearInFrontSpace(size_t space, size_t ind_inx, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> samplePush(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleByIndPos(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd);

		std::vector<std::vector<Real>> sampleByPush(Solution<>& sol, const std::vector<std::shared_ptr<Solution<>>>& his_sols, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleAmongGlobal(Population<Solution<>>& pop, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleAmongNeighs(size_t method, std::map<size_t, std::vector<size_t>>& pop2space, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleAmongNeighs(Population<Solution<>>& pop, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleAmongNeighs(Population<Solution<>>& pop, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleAmongSubspace(Population<Solution<>>& pop, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleAmongSubspace(Solution<>& ind, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleRandom(Solution<>& ind, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);

		bool spaceSubdivision(Problem* pro, Random* rnd);

		//void generateWeighVectors(size_t num_weigh);
		void selectNewPop();
		void ensembleSelect(size_t select_num, Problem* pro, Random* rnd);
		void SolutionSelection(Problem* pro, Random* rnd) override;
		void doubleClusterSelection(size_t select_num, Problem* pro, Random* rnd);
		void subspaceClusterSelection(size_t select_num, Problem* pro, Random* rnd);
		void objspaceClusterSelection(size_t select_num, Problem* pro, Random* rnd, bool add_his);
		void clusterSelection(size_t select_num, Problem* pro, Random* rnd, size_t m_num_neigh);

		void selectInHisRange(size_t select_num, Problem* pro, Random* rnd);

		bool indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol);

		void adaptiveSplit(Problem* pro);
		void adaptiveSplitFrontSpace(Problem* pro, Random* rnd);
		void adaptiveSplitSpace(Problem* pro, Random* rnd);
		void periodSplitSpace(Problem* pro, Random* rnd);
		//std::vector<std::pair<size_t,bool>> subspaceSeperable(size_t inx);//检测子空间是否具有可分性
		void splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd);
		//void splitSubspace(size_t inx,size_t num, int flag);
		void splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag);
		int findSplitDim(int inx, Problem* pro);

		void findClusterCenterSsp();
		void updateLinkSubspace(size_t inx, const std::vector<size_t>& front_spaces, Random* rnd);

		void testCoverage(Problem* pro);
		void updateNeighSpace();



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
		//std::unique_ptr<Subspace> m_exploit_box;
		std::unique_ptr<SpaceDirectedGraph> m_subspace_graph;
		bool m_update_cluster = false;
		size_t m_num_series = 10;//用于环境变化时间序列预测的长度
		std::vector<std::pair<Real, Real>> m_pop_range;//当前种群范围
		std::vector<std::pair<Real, Real>> m_pop_var_range;//当前种群范围
		std::vector<std::pair<Real, Real>> m_front_pop_range;//前沿种群目标范围
		//std::vector <std::shared_ptr<Solution<>>> m_history_front_sols;//store all non-diminated Solutions
		std::vector <std::shared_ptr<Solution<>>> m_archive;//store a number of front sols
		std::vector <std::shared_ptr<Population<Solution<>>>> m_gen_front_sols;//store front sols in each generation

		std::vector<Solution<>> m_subobj_opt_sol;//子目标最优的边界个体
		//std::vector<std::unique_ptr<ObjRegionInfo>> m_obj_region_info; //目标子空间信息
		size_t m_converge_age = 10;//设置收敛与否的阈值

		std::vector<size_t> m_over_space_index;//记录选择个体数不足的子空间
		std::vector<std::vector<size_t>> m_front_space_inx;//前沿子空间索引

		bool m_normalize = true;//目标空间归一化标记
		bool m_add_neighbor = false;//开发吸引域聚类时候含前排子空间邻域
		bool m_sample_in_basin = false;//基于子空间内采样还是基于吸引域采样
		bool m_subspace_select = false;//环境选择方式
		bool m_evolve_by_potential;//子代生成基于潜力采样还是传统方式
		bool m_mean_potential = true;//吸引域潜力为均值还是最值
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

		std::vector<std::map<size_t, Real>> m_all_front_nei_prob;//子空间连续性概率

		size_t m_total_push_times = 0;
		std::vector<size_t> m_push_results = { 0,0,0 };//dominant,dominated,nondominated

		size_t m_total_extend_times = 0;
		std::vector<size_t> m_extend_results = { 0,0,0 };

		std::vector<size_t> m_total_results = { 0,0,0 };//累积演化统计

		size_t m_push_behind_empty = 0;

		size_t m_queue_size = 20;

		Real m_test_rate;
		Real m_greedy_rate;
		Real m_disturbance_rate;
		Real m_neigh_interactive_rate;
		Real m_subspace_interactive_rate;

		size_t m_exploit_fail_times = 0;
		size_t m_exploit_succ_times = 0;

		bool m_exploit_flag = false;
		bool m_local_selection = false;
		std::vector<size_t> m_exploit_ind_index;
		std::vector<std::vector<std::pair<Real, Real>>> m_exploit_space_bounds;
		std::vector<std::vector<Real>> m_exploit_pos;

		//std::vector<std::pair<Real, Real>> m_reset_box;

		std::map<size_t, std::vector<size_t>> m_neigh_spaces;

		std::vector<size_t> m_ind_operator_times = { 0,0,0,0,0,0,0,0 };
		std::vector<std::vector<size_t>> m_operator_results;

		std::map<size_t, std::vector<size_t>> m_subspace_pop_index;

		std::vector<std::vector<size_t>> m_pop_index_clusters;
		std::vector<std::vector<size_t>> m_clusters_in_var;
		std::vector<std::vector<size_t>> m_pop_index;

		std::map<size_t, std::pair<size_t, std::queue<int>>> m_subspace_succ_rate;//space_idx,count,succ_rate

		std::map<size_t, std::pair<size_t, std::queue<int>>> m_subspace_interactive_result;//space_idx,count,succ_rate
	};

}

#endif  // !OFEC_SPMOEA1_H