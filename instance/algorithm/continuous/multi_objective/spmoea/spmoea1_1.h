/********* Begin Register Information **********
{
	"name": "SPMOEA1_1",
	"identifier": "SPMOEA1_1",
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
// 在spmoea1_16的基础上，根据子空间内的前沿解的线性度划分子空间

#ifndef OFEC_SPMOEA1_1_H
#define OFEC_SPMOEA1_1_H

#include "spmoea.h"

namespace ofec {

	class SPMOEA1_1 :public SPMOEA {
	public:

		void initialize_();
		void initPop(Problem* pro, Algorithm* alg, Random* rnd);
		void initiObjSpace(Problem* pro);
		void initiVarSpace(Problem* pro);
		void record() override;

		int evolve(Problem* pro, Algorithm* alg, Random* rnd);

		void calFrontBoundRatio(std::vector<Real>& front_bound_ratio);
		void PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro) override;
		void generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, const std::vector<size_t>& pop_resource, std::vector<int> type);

		bool spaceSubdivision(size_t space_inx, Problem* pro, Random* rnd);
		bool spaceSubdivision(Problem* pro, Random* rnd);

		bool indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol);

		int subspaceSeperable(size_t inx, Problem* pro) override;//检测子空间是否具有可分性
		void splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd);
		//void splitSubspace(size_t inx,size_t num, int flag);
		void splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag);
		int findSplitDim(int inx, Problem* pro);
		

		void findClusterCenterSsp();

		void updateFrontRegionLinkSpace(Problem* pro,Random* rnd) override;

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
		size_t m_divide_granularity=20;
		

	};

}

#endif  // !OFEC_SPMOEA1_1_H