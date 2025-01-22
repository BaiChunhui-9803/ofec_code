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
			size_t subspace_inx=0;                         //�ӿռ������
			std::list<size_t> m_sub_neighbors;             //�ӿռ������

			std::vector<Real> m_subObj_optima;             //�ӿռ���ÿ����Ŀ�������ֵ
			std::vector<size_t> m_subObj_optima_rank;      //�ӿռ���ÿ����Ŀ�������ֵ������ֵ
			std::vector<Real> m_subObj_improve;            //ÿ����Ŀ���������
			Real m_best_rank=INT16_MAX;                     //�ӿռ��ڵ�ǰ��õ�rankֵ
			Real m_sub_best_rank = INT16_MAX;               //�ӿռ�������Ⱥ�������ڵ����rankֵ
			size_t m_subspace_granularity=1;               //�ӿռ�����
			
			std::vector<size_t> m_linear_neigh_space_inx;  //���������ӿռ�����
			std::vector<std::shared_ptr<Solution<>>> m_history_inds;//�ӿռ��ڵ���ʷ��
			std::vector<std::shared_ptr<Solution<>>> m_curr_sols;      //�ӿռ��ڵ�ǰ�ĸ���
			std::vector<std::shared_ptr<Solution<>>> m_subspace_front_sol; //ÿ���ӿռ��ڵ�ǰ�ؽ�
			std::vector<size_t> m_subspace_front_inx;  //�ӿռ�ǰ�ؽ�����ʷ���е�����
			std::vector<std::shared_ptr<Solution<>>> m_represent_sol;           //�ӿռ��ڵĴ����
			std::vector<size_t> m_subspace_represent_inx;  //�ӿռ���������ʷ���е�����
			std::vector<std::shared_ptr<Solution<>>> m_front_sol_in_subspace; //�ӿռ�����ǰ���ϵĽ�
			std::vector<std::shared_ptr<Solution<>>> m_history_sub_obj_sol;//�ӿռ��ڵ���Ŀ����ʷ��ý�
			std::vector <std::shared_ptr<Population<Solution<>>>> m_gen_front_sols;//store front sols in each generation
			std::vector<int> idx_cluster;//����������
			size_t m_sub_freq = 0;   //�ӿռ����Ƶ��
			Real m_coverage_ratio = 0.; //�ӿռ�ĸ�����
			size_t m_linear_flag = 0;//0:nonlinear, 1:linear
			std::vector<size_t> m_link_subspaces;//����������

			std::vector<std::shared_ptr<Solution<>>> m_subspace_keep_pop;//�ӿռ��ڱ������ʱ��Ⱥ

			size_t m_represent_num = 3;


			Real m_sparse_value=0.; //ǰ���ӿռ��ϡ��ֵ
			Real m_cover_rate=1.; //�ӿռ�ĸ����ʣ����ø���֮����С����ֵ����

			Real m_potential = 0.;
			Real m_explore_potential = 0.;       //�ӿռ�̽��Ǳ��ֵ
			Real m_exploit_potential = 0.;        //�ӿռ俪��Ǳ��ֵ

			//int aff_bsn_atr = 0;

			//Real m_coff_sample = 1.;//����ϵ��
			Real m_coff_feedback = 1.;//����ϵ��
			
			size_t m_front_times = 0;//�ӿռ���ǰ�صĴ���
			//std::vector<size_t> m_rank_series;//�ӿռ�rankֵ�ı仯����
			size_t fail_times=0;//�ӿռ�����̽��ʧ�ܵĴ���
			bool improve_flag = true;    //�����Ƿ�����,ǰ��rankֵ�ı仯
			std::string search_flag = "explore";//����̽�����ǿ����ӿռ�
			std::string add_flag = "no";//�����ӿռ��Ƿ��оֲ�ǰ��

			bool m_update_flag = false;//�ӿռ���Ϣ���±��

			int m_feature_state = 1;//0��û�оֲ��ṹ��1���оֲ��ṹ
		};
		struct BasinInfo_MOP {
			size_t m_best_rank;//����������õ�����ֵ
			std::vector<std::shared_ptr<Solution<>>> m_current_indi;//�������ڵ�ǰ�ĸ���
			std::vector<size_t> m_local_front_solutions;//���������ھֲ�ǰ���ϵĽ�
			std::vector<size_t> m_front_solution;//����������ǰ���ϵĽ�
			int m_cluster_inx = -1;//�����������
			
			Real m_basin_potential;//�������̽��Ǳ��
			Real m_basin_explore_potential;//�������̽��Ǳ��
			Real m_basin_exploit_potential;//������Ŀ���Ǳ��
			size_t m_basin_freq;//�����������Ƶ��
			std::vector<size_t> m_subspace_set;//������������ӿռ�
			std::vector<std::shared_ptr<Solution<>>> m_history_inds;//�������ڵ���ʷ��
			Real m_coff_sample = 1.;//�������ڲ���ϵ��
			Real m_coff_feedback = 1.;//�������ڲ�������ϵ��
			size_t fail_times;//����������̽��ʧ�ܵĴ���
			std::string flag = "explore";
		};

		struct ObjRegionInfo {
			size_t subspace_inx = 0;                         //�ӿռ������
			std::list<size_t> m_obj_neighbors;             //�ӿռ������
			std::vector<Real> obj_optima;//�ӿռ���ÿ����Ŀ�������ֵ
			std::vector<std::shared_ptr<Solution<>>> m_curr_sols;//�ӿռ��ڵĸ��������,ģ������޸�Ϊָ����ָ��
			int m_obj_rank = std::numeric_limits<int>::max();//�ӿռ������ֵ
			//size_t m_duration = 0;//�����н�Ĵ���
			bool converged_region = false;//�Ƿ�Ϊ�������ӿռ�
			std::vector<size_t> ind_idx;//�ӿռ�����ĸ��������
		};

	private:
		std::unique_ptr<KDTree_MOP> m_split_tree_MOP;//��Ŀ�������ӿռ�
		std::vector<std::unique_ptr<SubspaceInfo_MOP>> m_subspace_info;
		std::vector<std::vector<std::vector<size_t>>> m_clusters_MOP;//���������������
		std::vector<std::vector<size_t>> m_front_clusters;//ǰ���ӿռ��������
		std::vector<std::vector<size_t>> m_cluster_centers;//���������ӿռ�����
		std::vector<std::unique_ptr<BasinInfo_MOP>> m_basin_info;
		//std::vector<std::unique_ptr<HLC>> m_hlc_subobj;//��Ŀ���HLC���Ƿ��б�Ҫ����Ϊ�Ѿ�������ÿ����Ŀ���������ӿռ�

		std::unique_ptr<KDTree_MOP> m_obj_tree_MOP;//��Ŀ�������ӿռ�
		std::vector<std::unique_ptr<ObjRegionInfo>> m_obj_space_info;
		int m_layer_num;
		

		std::unique_ptr<KDTree_MOP> m_selection_tree_MOP;//��Ŀ�������ӿռ�ѡ��

		Real m_total_volume;//�����ռ����

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
		void rankClustering();//����rankֵ����
		int getLayerNum() { return m_layer_num; };
		void setLayerNum(int v) { m_layer_num = v; }

		Real getSpaceVolume() { return m_total_volume; }
		void setSpaceVoume(Real v) { m_total_volume = v; }

		std::vector<std::vector<std::vector<size_t>>> clusterSubspace(const std::vector<size_t>& subspaces,size_t total_num_spaces,bool neighbor_flag);//����ĳһָ����ӿռ���о���
		std::vector<std::vector<size_t>> clusterFrontSpace(const std::vector<size_t>& frontspace, bool neighbor_flag);//ǰ���ӿռ����

		void updateSubspaceInfo(Population<Solution<>> &pop,Problem *pro, Random *rnd);//���ݵ�ǰ��ķֲ������ӿռ���Ϣ
		void updateSubSpaceInfo(std::vector<std::shared_ptr<Solution<>>>& pop, Problem *pro, Random *rnd);
		void updateSubspaceInfo(Solution<>& sol);//���ݵ�ǰ������ӿռ���Ϣ
		void updateBasinInfo(const std::vector<std::vector<size_t>>& explore_basin, const std::vector<std::vector<size_t>>& exploit_basin,bool b);//���ݵ�ǰ��ķֲ�������������Ϣ
		
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