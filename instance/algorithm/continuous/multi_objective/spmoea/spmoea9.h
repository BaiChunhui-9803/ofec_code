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

	typedef std::vector<std::vector<std::map<size_t, std::vector<size_t>>>> pop_attach;//Ŀ��������ռ������
	typedef std::vector<std::map<size_t, std::vector<size_t>>> space_attach;//��Ⱥ��ĳһ�ռ��������ǰ�źͷ�ǰ���ӿռ�
	using KDTree = nanoflann::KDTreeSpace<Real>;

	class SPMOEA9 :public Algorithm {
	public:
		struct ObjRegionInfo {
			std::vector<Real> obj_optima;//�ӿռ���ÿ����Ŀ�������ֵ
			std::vector<std::shared_ptr<Solution<>>> m_curr_sols;//�ӿռ��ڵĸ��������,ģ������޸�Ϊָ����ָ��
			int m_obj_rank = std::numeric_limits<int>::max();//�ӿռ������ֵ
			//size_t m_duration = 0;//�����н�Ĵ���
			bool converged_region = false;//�Ƿ�Ϊ�������ӿռ�
			std::vector<size_t> ind_idx;//�ӿռ�����ĸ��������
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

		space_attach spaceAttach(const Population<Solution<>>& pop, const KDTree& tree, bool b);//bѡ��Ŀ��ռ仹�������ռ�
		pop_attach popAttach(const Population<Solution<>>& pop);
		//std::vector<std::vector<std::vector<Real>>> indi_time_series(const std::vector<std::vector<Solution<>>>& p);
		//std::vector<Real>& optimaPredict(const std::vector<std::vector<Real>>& points);
		std::vector<size_t> selectIndi(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd);
		//std::vector<size_t> boundary_space(const std::vector<size_t>& index);//�ҳ�Ŀ��ռ�ı߽��ӿռ�
		std::vector<std::vector<size_t>> clusterExploitSpace(const KDTree& tree, const std::map<size_t, std::vector<size_t>>& frontspace);//��ǰ���ӿռ�����������
		std::vector<std::vector<size_t>> clusterExploreSpace(const KDTree& tree);//����ǰ���ӿռ䰴��̽������������
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
		//std::unique_ptr<KDTree> m_obj_tree; //Ŀ��ռ�ָ���
		std::unique_ptr<MO_HLC> m_mo_hlc;   //�����ռ�ָ���
		size_t m_pop_size;
		Real m_cr, m_mr, m_ceta, m_meta;
		size_t m_num_region_var = 0;
		size_t m_num_region_obj = 0;
		size_t m_num_series = 10;//���ڻ����仯ʱ������Ԥ��ĳ���
		std::vector<std::pair<Real, Real>> m_pop_range;//��ǰ��Ⱥ��Χ
		std::vector<std::pair<Real, Real>> m_pop_var_range;//��ǰ��Ⱥ��Χ
		std::vector<std::pair<Real, Real>> m_front_pop_range;//ǰ����ȺĿ�귶Χ
		std::vector <std::shared_ptr<Solution<>>> m_history_front_sols;//store all non-diminated Solutions
		std::vector <std::shared_ptr<Solution<>>> m_archive;//store a number of front sols

		std::vector <std::shared_ptr<Population<Solution<>>>> m_gen_front_sols;//store front sols in each generation

		std::vector<Solution<>> m_subobj_opt_sol;//��Ŀ�����ŵı߽����
		std::vector<std::unique_ptr<ObjRegionInfo>> m_obj_region_info; //Ŀ���ӿռ���Ϣ
		size_t m_converge_age = 10;//��������������ֵ
		Real m_R1,m_R2,m_R3;//����������ǰ��Ⱥ��Ŀ��ռ�ı�ֵ
		Real m_IGD;
		std::vector<size_t> m_over_space_index;//��¼ѡ�������������ӿռ�
		std::vector<std::vector<size_t>> m_front_space_inx;//ǰ���ӿռ�����
		
		size_t m_archive_num;   //�浵�Ĺ�ģ
		bool m_normalize = true;//Ŀ��ռ��һ�����
		bool m_add_neighbor = false;//�������������ʱ��ǰ���ӿռ�����
		bool m_sample_in_basin = false;//�����ӿռ��ڲ������ǻ������������
		bool m_subspace_select = false;//����ѡ��ʽ
		bool m_evolve_by_potential = true;//�Ӵ����ɻ���Ǳ���������Ǵ�ͳ��ʽ
		bool m_mean_potential = true;//������Ǳ��Ϊ��ֵ������ֵ
		std::vector<Real> m_update_tree_flag;
		std::vector<std::vector<Real>> m_num_e_e;
	};

}

#endif  // !OFEC_SPMOEA9_H
