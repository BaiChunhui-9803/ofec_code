/********* Begin Register Information **********
{
	"name": "SPMOEA3",
	"identifier": "SPMOEA3",
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
// Created: 16 Feb 2022 by Qingshan Tan (email:qingshan.t@cug.edu.cn)
// split space as pop evolve
// record information in each population
// ʶ���ӿռ��еľֲ��ṹ�����֣�����
// identify local structure and seperate them, then restart in each sub-box.
// k-d tree is constructed to record infos in each subspace for each sub-box.
// whether stop explore local structure or not according to the sampling uniformity in each sub-box

#ifndef OFEC_SPMOEA3_H
#define OFEC_SPMOEA3_H

#include "spmoea.h"

namespace ofec {

	class SPMOEA3 :public SPMOEA {
	public:
		void initialize_();
		void initPop(Problem* pro, Algorithm* alg, Random* rnd);
		void initiObjSpace(Problem* pro);
		void initiVarBox(std::vector<std::pair<Real, Real>>& box, Problem* pro);
		void record() override;

		int evolve(Problem* pro, Algorithm* alg, Random* rnd);

		void generatePop(std::vector<std::pair<Real, Real>>& bound, Problem* pro, Algorithm* alg, Random* rnd);
		void generateOffspring(Problem* pro, Random* rnd);
		void updatePop(size_t inx, std::vector<std::pair<Real, Real>>& bound, Problem* pro, Algorithm* alg, Random* rnd);

		void updateHistoryFrontSol(Population<Solution<>>& pop, Problem* pro);
		void updateArchive(size_t num, Problem* pro);

		int subspaceSeperable(size_t inx);//����ӿռ��Ƿ���пɷ���
		void splitSubspace(size_t inx, int flag);
		std::vector<std::pair<Real, Real>> distFromBound(size_t pop_inx);
		std::vector<Real> minDistFromBound(size_t pop_inx);

		bool updateCluster();
		void clusterSubspace();

		void setVarPartitionNum(const size_t v) { m_num_region_var = v; }
		void setObjPartitionNum(const size_t v) { m_num_region_obj = v; }
		size_t getVarRegionNum() const { return m_num_region_var; }
		size_t getObjRegionNum() const { return m_num_region_obj; }
		//MO_HLC& getMO_HLC(size_t inx) const { return *m_multi_hlc[inx]; }

#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;

	private:
		//std::vector<std::unique_ptr<MO_HLC>> m_multi_hlc;   //�ռ�ָ���
		size_t m_pop_size;
		Real m_cr, m_mr, m_ceta, m_meta;
		//std::unique_ptr<MultiPopulation<SPMOEA_pop>> m_multi_pops;
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
		size_t m_converge_age = 10;//��������������ֵ
		Real m_R1, m_R2, m_R3;//����������ǰ��Ⱥ��Ŀ��ռ�ı�ֵ
		Real m_IGD;
		std::vector<size_t> m_over_space_index;//��¼ѡ�������������ӿռ�
		std::vector<std::vector<size_t>> m_front_space_inx;//ǰ���ӿռ�����

		bool m_normalize = true;//Ŀ��ռ��һ�����

		size_t m_max_split_times = 3;
		size_t m_cur_split_times = 0;
		size_t m_split_num = 2;

		size_t m_sub_pop_size = 10;
		//�ӿռ䵥���Դ��
	};

}

#endif  // !OFEC_SPMOEA3_H

