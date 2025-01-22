/********* Begin Register Information **********
{
	"name": "SPMOEA5",
	"identifier": "SPMOEA5",
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
// Created: 25 Feb 2022 by Qingshan Tan (email:qingshan.t@cug.edu.cn)
// single-objective local optima dig
// record information in each population

#ifndef OFEC_SPMOEA5_H
#define OFEC_SPMOEA5_H

#include "spmoea.h"

namespace ofec {

	class SPMOEA5 :public SPMOEA {
	public:
		void initialize_();
		void initPop(Problem *pro, Algorithm *alg, Random *rnd);
		void initiObjSpace(Problem *pro);
		void initiVarSpace(Problem *pro);
		void record() override;

		int evolve(Problem *pro, Algorithm *alg, Random *rnd);

		void generatePop(std::vector<std::pair<Real, Real>>& bound, Problem *pro, Algorithm *alg, Random *rnd);
		void generateOffspring(Problem *pro, Random *rnd);

		void updateHistoryFrontSol(Population<Solution<>>& pop, Problem *pro);
		void updateArchive(size_t num, Problem *pro);

		std::vector<std::pair<size_t, bool>> subspaceSeperable(size_t inx);//检测子空间是否具有可分性
		void splitSubspace(size_t inx);
		std::vector<std::pair<Real, Real>> distFromBound(size_t pop_inx);

		bool updateCluster();
		void clusterSubspace();

		void setVarPartitionNum(const size_t v) { m_num_region_var = v; }
		void setObjPartitionNum(const size_t v) { m_num_region_obj = v; }
		size_t getVarRegionNum() const { return m_num_region_var; }
		size_t getObjRegionNum() const { return m_num_region_obj; }


#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;

	private:
		size_t m_pop_size;
		Real m_cr, m_mr, m_ceta, m_meta;
		size_t m_num_region_var = 0;
		size_t m_num_region_obj = 0;
		size_t m_num_series = 10;//用于环境变化时间序列预测的长度
		std::vector<std::pair<Real, Real>> m_pop_range;//当前种群范围
		std::vector<std::pair<Real, Real>> m_pop_var_range;//当前种群范围
		std::vector<std::pair<Real, Real>> m_front_pop_range;//前沿种群目标范围
		std::vector <std::shared_ptr<Solution<>>> m_history_front_sols;//store all non-diminated Solutions
		std::vector <std::shared_ptr<Solution<>>> m_archive;//store a number of front sols
		std::vector <std::shared_ptr<Population<Solution<>>>> m_gen_front_sols;//store front sols in each generation

		std::vector<Solution<>> m_subobj_opt_sol;//子目标最优的边界个体
		size_t m_converge_age = 10;//设置收敛与否的阈值
		Real m_R1, m_R2, m_R3;//三个衡量当前种群在目标空间的比值
		Real m_IGD;
		std::vector<size_t> m_over_space_index;//记录选择个体数不足的子空间
		std::vector<std::vector<size_t>> m_front_space_inx;//前沿子空间索引

		bool m_normalize = true;//目标空间归一化标记

		size_t m_max_split_times = 3;
		size_t m_split_num = 5;

		size_t m_sub_pop_size = 10;
		//子空间单调性打标
	};

}

#endif  // !OFEC_SPMOEA5_H

