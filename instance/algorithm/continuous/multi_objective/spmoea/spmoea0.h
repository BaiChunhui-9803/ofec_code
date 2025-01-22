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
		
		//std::vector<size_t> boundary_space(const std::vector<size_t>& index);//�ҳ�Ŀ��ռ�ı߽��ӿռ�
		std::vector<std::vector<size_t>> clusterExploitSpace(const KDTree& tree, const std::map<size_t, std::vector<size_t>>& frontspace);//��ǰ���ӿռ�����������
		std::vector<std::vector<size_t>> clusterExploreSpace(const KDTree& tree);//����ǰ���ӿռ䰴��̽������������
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

		//space_attach spaceAttach(const Population<Solution<>>& pop, const KDTree& tree, bool b);//bѡ��Ŀ��ռ仹�������ռ�
		//pop_attach popAttach(const Population<Solution<>>& pop);
		std::vector<size_t> selectIndi(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd);

#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;

	private:
		//std::unique_ptr<SPMOEA_pop> m_pop;
		//std::unique_ptr<KDTree> m_obj_tree; //Ŀ��ռ�ָ���
		//std::unique_ptr<MO_HLC> m_mo_hlc;   //�����ռ�ָ���
		//size_t m_pop_size;
		//Real m_cr, m_mr, m_ceta, m_meta;
		
		size_t m_num_series = 10;//���ڻ����仯ʱ������Ԥ��ĳ���
		//std::vector<std::pair<Real, Real>> m_pop_range;//��ǰ��Ⱥ��Χ
		std::vector<std::pair<Real, Real>> m_pop_var_range;//��ǰ��Ⱥ��Χ
		//std::vector<std::pair<Real, Real>> m_front_pop_range;//ǰ����ȺĿ�귶Χ
		std::vector <std::shared_ptr<Solution<>>> m_history_front_sols;//store all non-diminated Solutions
		std::vector <std::shared_ptr<Solution<>>> m_archive;//store a number of front sols

		std::vector <std::shared_ptr<Population<Solution<>>>> m_gen_front_sols;//store front sols in each generation

		std::vector<size_t> m_front_subspaces;//�����ӿռ���Ϣ֮ǰ��ǰ���ӿռ�

		std::vector<Solution<>> m_subobj_opt_sol;//��Ŀ�����ŵı߽����
		//std::vector<std::unique_ptr<ObjRegionInfo>> m_obj_region_info; //Ŀ���ӿռ���Ϣ
		size_t m_converge_age = 10;//��������������ֵ
		Real m_R1, m_R2, m_R3;//����������ǰ��Ⱥ��Ŀ��ռ�ı�ֵ
		Real m_IGD;
		std::vector<size_t> m_over_space_index;//��¼ѡ�������������ӿռ�
		std::vector<std::vector<size_t>> m_front_space_inx;//ǰ���ӿռ�����

		//size_t m_archive_num;   //�浵�Ĺ�ģ
		bool m_normalize = true;//Ŀ��ռ��һ�����
		bool m_add_neighbor = false;//�������������ʱ��ǰ���ӿռ�����
		bool m_sample_in_basin = false;//�����ӿռ��ڲ������ǻ������������
		bool m_subspace_select = false;//����ѡ��ʽ
		bool m_evolve_by_potential = true;//�Ӵ����ɻ���Ǳ���������Ǵ�ͳ��ʽ
		bool m_mean_potential = true;//������Ǳ��Ϊ��ֵ������ֵ
		std::vector<Real> m_update_tree_flag;
		std::vector<std::vector<Real>> m_num_e_e;

		size_t m_exploit_pop_size;
		size_t m_explore_pop_size;
	};

}

#endif  // !OFEC_SPMOEA0_H
