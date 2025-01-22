/********* Begin Register Information **********
{
	"name": "SPMOEA0",
	"identifier": "SPMOEA0",
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
// Created: 20 Feb 2023 by Qingshan Tan (email:qingshan.t@cug.edu.cn)
// Last modified: 
// partition space and find front subspace, calculate potential of subspaces, assign 
// resource for subspaces
// focus on exploration

#ifndef OFEC_SPMOEA0_H
#define OFEC_SPMOEA0_H

#include "spmoea.h"

namespace ofec {

	class SPMOEA0 :public SPMOEA {
	public:
		void initialize_() override;
		void initPop(Problem *pro, Algorithm *alg, Random *rnd);
		void generatePop(const std::vector<std::vector<size_t>>& space_clusters, Problem *pro, Algorithm *alg, Random *rnd);
		void generateGoodPop(const std::vector<size_t>& space_clusters, Problem *pro, Algorithm *alg, Random *rnd);
		void generateOffspring(Problem *pro, Random *rnd);
		void initiObjSpace(Problem *pro);
		std::vector<std::pair<size_t, size_t>> frontChanged();
		std::vector<std::vector<size_t>> spaceCluster();
		void record() override;

		std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> predictObjSpacePotential(Problem *pro);
		void calSpacePotential(std::vector<size_t>& spaces, Problem *pro);
		std::vector<size_t> assignBasinResource(const std::vector<std::vector<size_t>>& potential_inds);
		//std::vector<size_t> predictVarSpacePotential(Problem *pro);
		//void updateVarSpacePotential();
		//void assignBasinResource();

		std::vector<std::vector<size_t>> clusterSubspaces();
		bool isLocalParetoRegion(size_t inx,Problem *pro,Algorithm *alg,Random *rnd);
		Real shiftPlane(std::vector<std::vector<Real>>& inputdata, std::vector<Real> weigh, bool flag);
		
		void updateMultiPop(const std::vector<std::vector<size_t>> &clusters, const std::vector<size_t> &pre_front, Problem *pro, Algorithm *alg, Random *rnd);

		std::vector<size_t> getPreFrontSubspace() { return m_front_subspaces; }
		void setPreFrontSubspace(std::vector<size_t> ssp) { m_front_subspaces=ssp; }

		//bool setLinks(std::vector<std::vector<Real>> &v1, std::vector<std::vector<Real>>& v2);
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
		

		std::vector<size_t> getExploitSpace();
		std::vector<size_t> getExploreSpace();

		
		//std::vector<std::vector<std::vector<Real>>> indi_time_series(const std::vector<std::vector<Solution<>>>& p);
		//std::vector<Real>& optimaPredict(const std::vector<std::vector<Real>>& points);
		
		//std::vector<size_t> boundary_space(const std::vector<size_t>& index);//找出目标空间的边界子空间
		std::vector<std::vector<size_t>> clusterExploitSpace(const KDTree& tree, const std::map<size_t, std::vector<size_t>>& frontspace);//将前沿子空间进行邻域聚类
		std::vector<std::vector<size_t>> clusterExploreSpace(const KDTree& tree);//将非前沿子空间按照探索吸引力聚类
		//std::vector<size_t> calBasinPotential(const std::vector<std::vector<size_t>>& basin);

		//Real calSubspaceSpan(const std::vector<Real>& ideal_point, const std::vector<Real>& nadir_point, const std::vector<std::vector<Real>>& points);
		//std::tuple<std::vector<size_t>, bool> selectIndSubspace(const std::vector<Real>& ref_point, const std::vector<size_t>& ind_index, int num);

		
		void recordHisFront(Problem *pro);
		//void recordMetrics(Problem *pro);



		std::vector<Real> getObjUpdateFre() { return m_update_tree_flag; }
		std::vector<std::vector<Real>> getEERatio() { return m_num_e_e; }
		//MO_HLC& getMO_HLC() const { return *m_mo_hlc; }

		int evolve(Problem *pro, Algorithm *alg, Random *rnd);
		//void generateOffspring(const std::vector<size_t>& resource, const std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>>& potential_inds, Problem *pro, Random *rnd);
		std::vector<Real> localMutationSol(const std::vector<Real>& sol, const std::vector<std::pair<Real, Real>>& boundary, Problem *pro, Random *rnd);
		std::vector<Real> vectorMutationSol(const std::vector<Real>& sol1, const std::vector<Real>& sol2, Problem *pro, Random *rnd);

		void envirSelection(Problem *pro, Random *rnd);



		//void repairSol(std::vector<Real>& sol,std::vector<size_t> space_inx,Random *rnd);

		//space_attach spaceAttach(const Population<Solution<>>& pop, const KDTree& tree, bool b);//b选择目标空间还是搜索空间
		//pop_attach popAttach(const Population<Solution<>>& pop);
		std::vector<size_t> selectIndi(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd);

#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;

	private:
		//std::unique_ptr<SPMOEA_pop> m_pop;
		//std::unique_ptr<KDTree> m_obj_tree; //目标空间分割树
		//std::unique_ptr<MO_HLC> m_mo_hlc;   //搜索空间分割树
		//size_t m_pop_size;
		//Real m_cr, m_mr, m_ceta, m_meta;
		
		size_t m_num_series = 10;//用于环境变化时间序列预测的长度
		//std::vector<std::pair<Real, Real>> m_pop_range;//当前种群范围
		std::vector<std::pair<Real, Real>> m_pop_var_range;//当前种群范围
		//std::vector<std::pair<Real, Real>> m_front_pop_range;//前沿种群目标范围
		std::vector <std::shared_ptr<Solution<>>> m_history_front_sols;//store all non-diminated Solutions
		std::vector <std::shared_ptr<Solution<>>> m_archive;//store a number of front sols

		std::vector <std::shared_ptr<Population<Solution<>>>> m_gen_front_sols;//store front sols in each generation

		std::vector<size_t> m_front_subspaces;//更新子空间信息之前的前沿子空间

		std::vector<Solution<>> m_subobj_opt_sol;//子目标最优的边界个体
		//std::vector<std::unique_ptr<ObjRegionInfo>> m_obj_region_info; //目标子空间信息
		size_t m_converge_age = 10;//设置收敛与否的阈值
		Real m_R1, m_R2, m_R3;//三个衡量当前种群在目标空间的比值
		Real m_IGD;
		std::vector<size_t> m_over_space_index;//记录选择个体数不足的子空间
		std::vector<std::vector<size_t>> m_front_space_inx;//前沿子空间索引

		//size_t m_archive_num;   //存档的规模
		bool m_normalize = true;//目标空间归一化标记
		bool m_add_neighbor = false;//开发吸引域聚类时候含前排子空间邻域
		bool m_sample_in_basin = false;//基于子空间内采样还是基于吸引域采样
		bool m_subspace_select = false;//环境选择方式
		bool m_evolve_by_potential = true;//子代生成基于潜力采样还是传统方式
		bool m_mean_potential = true;//吸引域潜力为均值还是最值
		std::vector<Real> m_update_tree_flag;
		std::vector<std::vector<Real>> m_num_e_e;

		size_t m_exploit_pop_size;
		size_t m_explore_pop_size;
	};

}

#endif  // !OFEC_SPMOEA0_H
