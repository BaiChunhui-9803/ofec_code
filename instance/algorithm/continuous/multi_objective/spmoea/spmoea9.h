/********* Begin Register Information **********
{
	"name": "SPMOEA9",
	"identifier": "SPMOEA9",
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
// Created: 28 Feb 2022 by Qingshan Tan (email:qingshan.t@cug.edu.cn)
// Last modified: 

#ifndef OFEC_SPMOEA9_H
#define OFEC_SPMOEA9_H

//#include "../../../template/selection/multi_objective/nsgaii.h"
#include "../../../../../core/algorithm/algorithm.h"
#include "../../../template/classic/ga/sbx_pop.h"
#include "../../../template/classic/de/population.h"
#include "../../../../../utility/kd-tree/kdtree_space.h"
#include "../mo_hlc.h"
#include "../../../../record/multi_objective/rcr_vec_real_moea.h"
#include "../../../../record/rcr_vec_real.h"
#include<tuple>

namespace ofec {

	class SPMOEA9_pop :public PopSBX<>{
	public:
		explicit SPMOEA9_pop(size_t size_pop, Problem *pro);
		void initialize(Problem *pro, Random *rnd) override;
		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
		void envirSelection(Problem *pro, Random *rnd, bool b);
		void NDSort(Population<Solution<>>& pop);
		Population<Solution<>>& getOffspring() { return m_offspring; }
		Population<Solution<>>& getParent() { return m_parents; }
		
	private:
		Population<Solution<>> m_offspring;  // 2 size of parent population
		Population<Solution<>> m_parents;   // parent pop 
	};

	typedef std::vector<std::vector<std::map<size_t, std::vector<size_t>>>> pop_attach;//目标和搜索空间的依附
	typedef std::vector<std::map<size_t, std::vector<size_t>>> space_attach;//种群在某一空间的依附，前排和非前排子空间
	using KDTree = nanoflann::KDTreeSpace<Real>;

	class SPMOEA9 :public Algorithm {
	public:
		struct ObjRegionInfo {
			std::vector<Real> obj_optima;//子空间内每个子目标的最优值
			std::vector<std::shared_ptr<Solution<>>> m_curr_sols;//子空间内的个体的索引,模板参数修改为指向解的指针
			int m_obj_rank = std::numeric_limits<int>::max();//子空间的排序值
			//size_t m_duration = 0;//持续有解的代数
			bool converged_region = false;//是否为收敛的子空间
			std::vector<size_t> ind_idx;//子空间包含的个体的索引
		};

		//SPMOEA()
		void initialize_() override;
		void initPop();
		void initiObjSpace(Problem *pro);
		void initiVarSpace(Problem *pro);
		void record() override;

		std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> predictObjSpacePotential(Problem *pro);
		void calSpacePotential(std::vector<size_t>& spaces,Problem *pro);
		std::vector<size_t> assignBasinResource(const std::vector<std::vector<size_t>>& potential_inds);
		//std::vector<size_t> predictVarSpacePotential(Problem *pro);
		//void updateVarSpacePotential();
		//void assignBasinResource();
		void updateObjSpaceInfo(const Population<Solution<>>& pop, Problem *pro, bool changed);
		bool updateObjRange(const Population<Solution<>>& pop, Problem *pro);
		bool updateHisFrontObjRange();
		void updateObjTree();
		void updateSolSpaceInfo(const Population<Solution<>>& pop, Problem *pro, bool b);
		void updateSubObjOpt(const Population<Solution<>>& pop);
		void updateSubObjRank(Problem *pro);
		void updateHistorySol(Problem *pro);
		void updateArchive(size_t num, Problem *pro);
		void updatePopAge();
		////void assignPop();
		//void updatePotential();
		bool popConverged();
		bool ifApproxiConverge(Random *rnd);
		void repairSol(std::vector<Real>& sol,Problem *pro, Random *rnd);

		std::vector<size_t> getExploitSpace();
		std::vector<size_t> getExploreSpace();

		space_attach spaceAttach(const Population<Solution<>>& pop, const KDTree& tree, bool b);//b选择目标空间还是搜索空间
		pop_attach popAttach(const Population<Solution<>>& pop);
		//std::vector<std::vector<std::vector<Real>>> indi_time_series(const std::vector<std::vector<Solution<>>>& p);
		//std::vector<Real>& optimaPredict(const std::vector<std::vector<Real>>& points);
		std::vector<size_t> selectIndi(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd);
		//std::vector<size_t> boundary_space(const std::vector<size_t>& index);//找出目标空间的边界子空间
		std::vector<std::vector<size_t>> clusterExploitSpace(const KDTree& tree, const std::map<size_t, std::vector<size_t>>& frontspace);//将前沿子空间进行邻域聚类
		std::vector<std::vector<size_t>> clusterExploreSpace(const KDTree& tree);//将非前沿子空间按照探索吸引力聚类
		//std::vector<size_t> calBasinPotential(const std::vector<std::vector<size_t>>& basin);
		
		//Real calSubspaceSpan(const std::vector<Real>& ideal_point, const std::vector<Real>& nadir_point, const std::vector<std::vector<Real>>& points);
		//std::tuple<std::vector<size_t>, bool> selectIndSubspace(const std::vector<Real>& ref_point, const std::vector<size_t>& ind_index, int num);

		void setVarPartitionNum(const size_t v) { m_num_region_var = v; }
		void setObjPartitionNum(const size_t v) { m_num_region_obj = v; }
		size_t getVarRegionNum() const { return m_num_region_var; }
		size_t getObjRegionNum() const { return m_num_region_obj; }
		std::vector<std::pair<Real, Real>> getPopRange() { return m_pop_range; }
		std::vector<std::pair<Real, Real>> getFrontPopRange() { return m_front_pop_range; }
		size_t numObjSubspace() { return m_obj_region_info.size(); }
		//const KDTree& getObjspaceTree() { return *m_obj_tree; }
		const ObjRegionInfo& getObjRegionInfo(size_t idx) { return *m_obj_region_info[idx]; }
		const SPMOEA9_pop& getPop() { return *m_pop; }
		void recordHisFront(Problem *pro);
		void recordMetrics(Problem *pro);
		


		std::vector<Real> getObjUpdateFre(){return m_update_tree_flag;}
		std::vector<std::vector<Real>> getEERatio() { return m_num_e_e; }
		MO_HLC& getMO_HLC() const { return *m_mo_hlc; }
		
		int evolve(Problem *pro, Algorithm *alg, Random *rnd);
		void generateOffspring(const std::vector<size_t>& resource,const std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>>& potential_inds,Problem *pro, Random *rnd);
		std::vector<Real> localMutationSol(const std::vector<Real>& sol, const std::vector<std::pair<Real, Real>>& boundary,Problem *pro, Random *rnd);
		std::vector<Real> vectorMutationSol(const std::vector<Real>& sol1, const std::vector<Real>& sol2,Problem *pro, Random *rnd);

		void envirSelection(Problem *pro, Random *rnd);

#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;
		
	private:
		std::unique_ptr<SPMOEA9_pop> m_pop;
		//std::unique_ptr<KDTree> m_obj_tree; //目标空间分割树
		std::unique_ptr<MO_HLC> m_mo_hlc;   //搜索空间分割树
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
		std::vector<std::unique_ptr<ObjRegionInfo>> m_obj_region_info; //目标子空间信息
		size_t m_converge_age = 10;//设置收敛与否的阈值
		Real m_R1,m_R2,m_R3;//三个衡量当前种群在目标空间的比值
		Real m_IGD;
		std::vector<size_t> m_over_space_index;//记录选择个体数不足的子空间
		std::vector<std::vector<size_t>> m_front_space_inx;//前沿子空间索引
		
		size_t m_archive_num;   //存档的规模
		bool m_normalize = true;//目标空间归一化标记
		bool m_add_neighbor = false;//开发吸引域聚类时候含前排子空间邻域
		bool m_sample_in_basin = false;//基于子空间内采样还是基于吸引域采样
		bool m_subspace_select = false;//环境选择方式
		bool m_evolve_by_potential = true;//子代生成基于潜力采样还是传统方式
		bool m_mean_potential = true;//吸引域潜力为均值还是最值
		std::vector<Real> m_update_tree_flag;
		std::vector<std::vector<Real>> m_num_e_e;
	};

}

#endif  // !OFEC_SPMOEA9_H
