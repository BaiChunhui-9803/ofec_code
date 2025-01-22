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
// learning from historical data in the search space for MOPs
// 
//  Created: 28 Feb 2022 by Qingshan Tan (email:qingshan.t@cug.edu.cn)
// Last modified: 

#ifndef OFEC_MO_HLC_H
#define OFEC_MO_HLC_H

//#include "../../template/framework/spae/hlc.h"
#include "../../../../utility/kd-tree/kdtree_space.h"
#include "../../../../core/algorithm/population.h"
#include "../../../../core/problem/solution.h"
//#include "solution_mop.h"
//#include <memory>
//#include <list>

namespace ofec {
	
	class MO_HLC {
	public:
		using KDTree_MOP = nanoflann::KDTreeSpace<Real>;

		struct SubspaceInfo_MOP {
			size_t subspace_inx=0;                         //子空间的索引
			std::list<size_t> m_sub_neighbors;             //子空间的邻域

			std::vector<Real> m_subObj_optima;             //子空间内每个子目标的最优值
			std::vector<size_t> m_subObj_optima_rank;      //子空间内每个子目标的最优值的排序值
			std::vector<Real> m_subObj_improve;            //每个子目标的提升量
			Real m_best_rank=INT16_MAX;                     //子空间内当前最好的rank值
			Real m_sub_best_rank = INT16_MAX;               //子空间在子种群吸引域内的最好rank值
			size_t m_subspace_granularity=1;               //子空间粒度
			
			std::vector<size_t> m_linear_neigh_space_inx;  //线性邻域子空间索引
			std::vector<std::shared_ptr<Solution<>>> m_history_inds;//子空间内的历史解
			std::vector<std::shared_ptr<Solution<>>> m_curr_sols;      //子空间内当前的个体
			std::vector<std::shared_ptr<Solution<>>> m_subspace_front_sol; //每个子空间内的前沿解
			std::vector<size_t> m_subspace_front_inx;  //子空间前沿解在历史解中的索引
			std::vector<std::shared_ptr<Solution<>>> m_represent_sol;           //子空间内的代表解
			std::vector<size_t> m_subspace_represent_inx;  //子空间代表解在历史解中的索引
			std::vector<std::shared_ptr<Solution<>>> m_front_sol_in_subspace; //子空间内在前沿上的解
			std::vector<std::shared_ptr<Solution<>>> m_history_sub_obj_sol;//子空间内的子目标历史最好解
			std::vector <std::shared_ptr<Population<Solution<>>>> m_gen_front_sols;//store front sols in each generation
			std::vector<int> idx_cluster;//所属吸引域
			size_t m_sub_freq = 0;   //子空间采样频率
			Real m_coverage_ratio = 0.; //子空间的覆盖率
			size_t m_linear_flag = 0;//0:nonlinear, 1:linear
			std::vector<size_t> m_link_subspaces;//所属吸引域

			std::vector<std::shared_ptr<Solution<>>> m_subspace_keep_pop;//子空间内保存的临时种群

			size_t m_represent_num = 3;


			Real m_sparse_value=0.; //前沿子空间的稀疏值
			Real m_cover_rate=1.; //子空间的覆盖率，采用个体之间最小距离值计算

			Real m_potential = 0.;
			Real m_explore_potential = 0.;       //子空间探索潜力值
			Real m_exploit_potential = 0.;        //子空间开发潜力值

			//int aff_bsn_atr = 0;

			//Real m_coff_sample = 1.;//采样系数
			Real m_coff_feedback = 1.;//反馈系数
			
			size_t m_front_times = 0;//子空间在前沿的次数
			//std::vector<size_t> m_rank_series;//子空间rank值的变化序列
			size_t fail_times=0;//子空间连续探索失败的次数
			bool improve_flag = true;    //采样是否提升,前后rank值的变化
			std::string search_flag = "explore";//属于探索还是开发子空间
			std::string add_flag = "no";//表征子空间是否有局部前沿

			bool m_update_flag = false;//子空间信息更新标记

			int m_feature_state = 1;//0：没有局部结构；1：有局部结构
		};
		struct BasinInfo_MOP {
			size_t m_best_rank;//吸引域内最好的排序值
			std::vector<std::shared_ptr<Solution<>>> m_current_indi;//吸引域内当前的个体
			std::vector<size_t> m_local_front_solutions;//吸引域内在局部前沿上的解
			std::vector<size_t> m_front_solution;//吸引域内在前沿上的解
			int m_cluster_inx = -1;//吸引域的索引
			
			Real m_basin_potential;//吸引域的探索潜力
			Real m_basin_explore_potential;//吸引域的探索潜力
			Real m_basin_exploit_potential;//吸引域的开发潜力
			size_t m_basin_freq;//吸引域的搜索频率
			std::vector<size_t> m_subspace_set;//吸引域包含的子空间
			std::vector<std::shared_ptr<Solution<>>> m_history_inds;//吸引域内的历史解
			Real m_coff_sample = 1.;//吸引域内采样系数
			Real m_coff_feedback = 1.;//吸引域内采样反馈系数
			size_t fail_times;//吸引域连续探索失败的次数
			std::string flag = "explore";
		};

		struct ObjRegionInfo {
			size_t subspace_inx = 0;                         //子空间的索引
			std::list<size_t> m_obj_neighbors;             //子空间的邻域
			std::vector<Real> obj_optima;//子空间内每个子目标的最优值
			std::vector<std::shared_ptr<Solution<>>> m_curr_sols;//子空间内的个体的索引,模板参数修改为指向解的指针
			int m_obj_rank = std::numeric_limits<int>::max();//子空间的排序值
			//size_t m_duration = 0;//持续有解的代数
			bool converged_region = false;//是否为收敛的子空间
			std::vector<size_t> ind_idx;//子空间包含的个体的索引
		};

	private:
		std::unique_ptr<KDTree_MOP> m_split_tree_MOP;//多目标搜索子空间
		std::vector<std::unique_ptr<SubspaceInfo_MOP>> m_subspace_info;
		std::vector<std::vector<std::vector<size_t>>> m_clusters_MOP;//聚类的吸引域索引
		std::vector<std::vector<size_t>> m_front_clusters;//前沿子空间聚类索引
		std::vector<std::vector<size_t>> m_cluster_centers;//聚类中心子空间索引
		std::vector<std::unique_ptr<BasinInfo_MOP>> m_basin_info;
		//std::vector<std::unique_ptr<HLC>> m_hlc_subobj;//子目标的HLC，是否有必要，因为已经保留了每个子目标最优在子空间

		std::unique_ptr<KDTree_MOP> m_obj_tree_MOP;//多目标搜索子空间
		std::vector<std::unique_ptr<ObjRegionInfo>> m_obj_space_info;
		int m_layer_num;
		

		std::unique_ptr<KDTree_MOP> m_selection_tree_MOP;//多目标搜索子空间选择

		Real m_total_volume;//搜索空间体积

	public:
		/*MO_HLC() {
			std::cout << "test" << std::endl;
		}*/
		
		MO_HLC(size_t dim_var, size_t num_obj);
	
		void initialVarSpace(std::vector<std::pair<Real, Real>>& boundary, int subspace_num);
		void initialVarSelectionSpace(std::vector<std::pair<Real, Real>>& boundary, int subspace_num);
		void initialObjSpace(std::vector<std::pair<Real, Real>>& boundary, int subspace_num);
		void initialBasin(size_t num_basin);
		void getSampling(const Solution<>& sample_point);

		void findClusterCenterSsp();
		void rankClustering();//基于rank值聚类
		int getLayerNum() { return m_layer_num; };
		void setLayerNum(int v) { m_layer_num = v; }

		Real getSpaceVolume() { return m_total_volume; }
		void setSpaceVoume(Real v) { m_total_volume = v; }

		std::vector<std::vector<std::vector<size_t>>> clusterSubspace(const std::vector<size_t>& subspaces,size_t total_num_spaces,bool neighbor_flag);//基于某一指标对子空间进行聚类
		std::vector<std::vector<size_t>> clusterFrontSpace(const std::vector<size_t>& frontspace, bool neighbor_flag);//前沿子空间聚类

		void updateSubspaceInfo(Population<Solution<>> &pop,Problem *pro, Random *rnd);//根据当前解的分布更新子空间信息
		void updateSubSpaceInfo(std::vector<std::shared_ptr<Solution<>>>& pop, Problem *pro, Random *rnd);
		void updateSubspaceInfo(Solution<>& sol);//根据当前解更新子空间信息
		void updateBasinInfo(const std::vector<std::vector<size_t>>& explore_basin, const std::vector<std::vector<size_t>>& exploit_basin,bool b);//根据当前解的分布更新吸引域信息
		
		void calBasinPotential(){}

		void converged(const Solution<>& center);

		//Real calSubspacePotential(size_t idx);
		void rankSubspaceByPotential(const std::vector<size_t>& subspaces);
		//Real calExploitSpacePotential(size_t idx,std::vector<size_t>& spaces);
		//Real calExploreSpacePotential(size_t idx);
		
		void spaceDivide(int idx, size_t sub_num);
		void spaceSplit(int idx, int dim, Real pos);

		void subspaceMerge(std::vector<int> idx){}

		//void update_subspace_potential();
		//void update_basin_potential();
		//std::vector<std::vector<size_t>>& get_neighbors_MOP() { return m_neighbors; }

		KDTree_MOP& getObjspaceTree() { return *m_obj_tree_MOP; }
		void updateObjTree() { m_obj_tree_MOP.reset(new KDTree_MOP()); }
		void updateVarSelectionTree() { m_selection_tree_MOP.reset(new KDTree_MOP()); }
		KDTree_MOP& getVarSelectionTree() { return *m_selection_tree_MOP; }

		size_t numObjspace() { return m_obj_space_info.size(); }
		const ObjRegionInfo& getObjRegionInfo(size_t idx) { return *m_obj_space_info[idx]; }

		std::vector<std::vector<size_t>> getCluster(size_t i) { return m_clusters_MOP[i]; }
		std::vector<std::vector<size_t>> getFrontClusters() { return m_front_clusters; }
		void setFrontClusters(const std::vector<std::vector<size_t>>& front_clusters) { m_front_clusters = front_clusters; }
		void setSpaceClusters(const std::vector<std::vector<std::vector<size_t>>>& space_clusters) { m_clusters_MOP = space_clusters; }
		size_t numSubspace() { return m_subspace_info.size(); }
		size_t numBasin() { return m_clusters_MOP.size(); }
		KDTree_MOP& subspaceTree() { return *m_split_tree_MOP; }
		SubspaceInfo_MOP& getSubspaceInfo(size_t idx) { return *m_subspace_info.at(idx); }
		void deleteSubspaceInfo(size_t inx) { m_subspace_info.erase(m_subspace_info.begin()+inx); }

		BasinInfo_MOP& getBasinInfo(size_t idx) { return *m_basin_info.at(idx); }
		const Real getSpaceBestRank(size_t inx)const { return m_subspace_info.at(inx)->m_best_rank; }
		void setSpaceBestRank(size_t inx, Real val) { m_subspace_info.at(inx)->m_best_rank = val; }
		const size_t getBasinBestRank(size_t inx)const { return m_basin_info.at(inx)->m_best_rank; }
		void setBasinBestRank(size_t inx, size_t val) { m_basin_info.at(inx)->m_best_rank = val; }
		std::vector<std::vector<std::vector<size_t>>>& getClusters() { return m_clusters_MOP; }
		const std::vector<Real> getSubspacePotential(size_t inx) const 
		{ 
			std::vector<Real> potential;
			potential.push_back(m_subspace_info.at(inx)->m_explore_potential);
			potential.push_back(m_subspace_info.at(inx)->m_exploit_potential);
			return potential;
		}
		void setSubspacePotential(size_t inx, std::vector<Real> val) 
		{ 
			m_subspace_info.at(inx)->m_explore_potential = val[0];
			m_subspace_info.at(inx)->m_exploit_potential = val[1];
		}
		const std::vector<Real> getBasinPotential(size_t inx) const
		{ 
			std::vector<Real> potential;
			potential.push_back(m_basin_info.at(inx)->m_basin_explore_potential);
			potential.push_back(m_basin_info.at(inx)->m_basin_exploit_potential);
			return potential;
		}
		void setBasinPotential(size_t inx, std::vector<Real> val)
		{ 
			m_basin_info.at(inx)->m_basin_explore_potential = val[0];
			m_basin_info.at(inx)->m_basin_exploit_potential = val[1];
		}
	};
}

#endif