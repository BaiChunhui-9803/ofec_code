/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Xiaofang Wu
* Email: changhe.lw@google.com Or wuxiaofang@cug.edu.cn
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 11 July 2019
// Last modified: 22 Sep 2020 by Xiaofang Wu (email:email:wuxiaofang@cug.edu.cn)

#ifndef OFEC_KDT_MOEA_H
#define OFEC_KDT_MOEA_H

#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/problem/solution.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../utility/kd-tree/kdtree_space.h"
#include "../../../template/selection/multi_objective/nsgaii.h"
#include "../../../../../utility/nonDominated_sorting/filter_sort.h"
#include "../../../../../utility/nonDominated_sorting/fast_sort.h"

namespace ofec {
	template<typename Population, typename Solution>
	class KDT_MOEA {
	public:
		using KDTree = nanoflann::KDTreeSpace<Real>;
		KDT_MOEA(const ParameterMap& v) : m_num_partition(v.get<int>("numPartition")), m_region(v.get<int>("numPartition")), kdtree(std::make_unique<KDTree>()) { }
		void set_kdtree(std::vector<std::pair<Real, Real>>& boundary,Problem *pro);
		void updateRegionLabel(std::vector<Solution>& vec,Problem *pro);
		void sort_region(const Population& pop, std::vector<int>& rank_dds_sa, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind, int pops, int pops_F1,Problem *pro);	//non-dominance sort of dds and sa
		void sort_region_predict(const std::vector<std::vector<Real>>& prediction, const Population& pop, std::vector<int>& rank_dds_sa, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind, int pops, int pops_F1,Problem *pro);	//non-dominance sort of dds and sa

		void sort_region_new(const std::vector<std::vector<Real>>& ref_point, std::vector<int>& rank_norm, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind,Problem *pro);

		void sort_region_norm_angle(const std::vector<std::vector<Real>>& ref_point, const Population& pop, std::vector<int>& rank_norm_angle, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind, int pops,Problem *pro);	//non-dominance sort of norm and angle (or distance)

		void sort_dd(std::vector<Solution>& offspring);
		void sort_dd_unify(std::vector<Solution>& offspring,Problem *pro);
		Real dd(Solution& ind1, Solution& ind2,Problem *pro);

		void sort_1_k_dominance(std::vector<Solution>& offspring, int m_k,Problem *pro);
		int dominance_1_k(const std::vector<std::vector<Real>*>& data, std::vector<int>& rank, const std::vector<OptimizeMode>& opt_mode, int m_k);
		Dominance objective_compare_1_k(const std::vector<Real>& a, const std::vector<Real>& b, const std::vector<OptimizeMode>& mode, int m_k);


		void cal_dds(const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind,Problem *pro);	//calculate dominace degree of subspace(dds)
		void get_sub(std::vector<std::pair<Real, Real>>& sub, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind, int k);
		Real cal_D_ij(const std::vector<std::pair<Real, Real>>& sub1, const std::vector<std::pair<Real, Real>>& sub2,Problem *pro);	//calculate dominace ratio based on objective difference(i.e. "check_dominance_ratio_obj_difference" in KDTree_MOEA_NSGAII.h)
		Real cal_sa(const Population& pop, const std::vector<Solution>& offspring, const std::vector<int>& region_to_sort_ind_k, int k, int pops,Problem *pro); //calculate smallest angle(sa)
		Real cal_cdi(const Population& pop, const std::vector<Solution>& offspring, const std::vector<int>& region_to_sort_ind_k, int k, int pops,Problem *pro); //calculate smallest angle(sa)/norm
		Real cal_cdi_predict(const std::vector<std::vector<Real>>& prediction, const Population& pop, const std::vector<Solution>& offspring, const std::vector<int>& region_to_sort_ind_k, int k, int pops, std::vector<std::vector<Real>>& observe_cdi_k,Problem *pro); //calculate smallest angle(sa)/norm based on prediction points

		bool find_reference_point(const std::vector<std::vector<Real>>& prediction, const std::vector<Real>& s_k, std::vector<Real>& z_star,Problem *pro);//find reference point (z_star) for subspace s_k 

		void update_min_max(const Population& pop,Problem *pro);
		void update_min_max(const std::vector<Solution>& offspring,Problem *pro);

		void update_min_max_predict(const std::vector<std::vector<Real>>& prediction,Problem *pro); //update m_min_max_predict

		KDTree& get_kdt() { return *kdtree; }
		std::vector<std::pair<Real, Real>>& min_max() { return m_min_max; }
		int size_region() const { return m_region_sort.size(); }
		std::vector<std::pair<int, Real>>& region() { return m_region_sort; }
		std::vector<std::pair<int, Real>>& region_dds() { return m_region_sort1; }
		std::vector<std::pair<int, Real>>& region_sa() { return m_region_sort2; }
		std::vector<std::pair<int, Real>>& region_norm() { return m_region_sort1; }
		std::vector<std::pair<int, Real>>& region_angle() { return m_region_sort2; }
		std::vector<std::pair<int, Real>>& region_rank() { return m_region_sort3; }
		std::vector<std::pair<int, Real>>& region_norm_only() { return m_region_sort_norm; }
		int get_num_partition() const { return m_num_partition; }
		void set_F1(int F1) { m_count_F1 = F1; }
		int get_F1() { return m_count_F1; }

	protected:
		std::unique_ptr<KDTree> kdtree;
		int m_num_partition;

		int m_count_F1 = 0; //record size of Solutions of rank = 0

		std::vector<std::vector<Real>> m_D;

		std::vector<std::pair<int, bool>> m_region;	//the first is subregion labeling,the second is used to mark whether there is an Solution or not.
		std::vector<std::pair<int, Real>> m_region_sort; //the first is subregion labeling,the second is rank value.
		std::vector<std::pair<int, Real>> m_region_sort1;//record dds or norm
		std::vector<std::pair<int, Real>> m_region_sort2;//record sa or angle
		std::vector<std::pair<int, Real>> m_region_sort3;//record rank after non-dominance sort (dds and sa) or (norm and angle)
		std::vector<std::pair<int, Real>> m_region_sort_norm;//record norm between z_star and x_k

		std::vector<std::pair<Real, Real>> m_min_max;

		std::vector<std::pair<Real, Real>> m_min_max_predict; //update by offspring and prediction point

	public:
		bool m_observe = true;//observe dds and sa in the first selection
		std::vector<std::vector<std::vector<Real>>> m_observe_cdi;
		std::vector<std::vector<std::vector<Real>>> m_observe_norm;
		std::vector<std::vector<std::vector<Real>>> m_observe_angle;
		int m_size_sort_first = 0;

		std::vector<std::vector<std::vector<Real>>> m_observe_angle_sub_inner;//observe angle in subspace inner

	};

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::set_kdtree(std::vector<std::pair<Real, Real>>& boundary,Problem *pro) {
		int M = CAST_CONOP(pro)->numberObjectives();
		kdtree.reset(new KDTree(std::vector<Real>(m_num_partition, 1.0 / m_num_partition), boundary));
		//kdtree->setInitBox(boundary);
		//kdtree->inputRatioData(std::vector<Real>(m_num_partition, 1.0 / m_num_partition));
		kdtree->buildIndex();

		/*output bound of kd-tree*/
		/*for (int i = 0; i < m_num_partition; ++i) {
			for (int j = 0; j < M; ++j) {
				std::cout << i << ": " << "(" << kdtree->get_box(i)[j].first << "," << kdtree->get_box(i)[j].second << ")" << "\t";
			}
			std::cout << std::endl;
		}*/
	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::updateRegionLabel(std::vector<Solution>& vec,Problem *pro) {
		if (m_region_sort.size() != 0) {
			std::vector<std::pair<int, Real>>().swap(m_region_sort);		//free m_region_sort
			std::vector<std::pair<int, Real>>().swap(m_region_sort1);		//free m_region_sort1
			std::vector<std::pair<int, Real>>().swap(m_region_sort2);		//free m_region_sort2
			std::vector<std::pair<int, Real>>().swap(m_region_sort3);		//free m_region_sort3
			std::vector<std::pair<int, Real>>().swap(m_region_sort_norm);	//free m_region_sort_norm
		}

		int M = CAST_CONOP(pro)->numberObjectives();
		for (int i = 0; i < m_num_partition; ++i) {
			m_region[i].first = i;
			m_region[i].second = false;
		}
		int k;
		for (int i = 0; i < vec.size(); ++i) {
			if (vec[i].getType() != -1)
				continue;
			int flag = true;
			k = kdtree->getRegionIdx(vec[i].objective());
			for (int s = 0; s < m_region_sort.size(); ++s) {
				if (k == m_region_sort[s].first) {
					flag = false;
					break;
				}
			}
			if (flag) {
				m_region_sort.push_back(std::pair<int, Real>(k, -1));	//make_pair<int,int>(k,-1)	
				m_region[k].second = true;
			}

		}
		m_region_sort1 = m_region_sort; //record dds
		m_region_sort2 = m_region_sort; //record sa
		m_region_sort3 = m_region_sort; //record rank after non-dominance sort
		m_region_sort_norm = m_region_sort; //record norm between z_star and x_k

		m_size_sort_first = m_region_sort.size();
	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::sort_region(const Population& pop, std::vector<int>& rank_dds_sa, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind, int pops, int pops_F1,Problem *pro) {
		cal_dds(offspring, region_to_sort_ind,pro);	//calculate dds
		int size_region = m_region_sort.size();
		if (rank_dds_sa.size() != size_region)
			rank_dds_sa.resize(size_region);
		//non-dominance sort of dds and sa
		std::vector<std::vector<Real>> data_dds_sa(size_region, std::vector<Real>(2));
		std::vector<std::vector<Real>*> data_dds_sa1(size_region);
		for (int k = 0; k < size_region; ++k) {
			data_dds_sa[k][0] = m_region_sort[k].second;
			//data_dds_sa[k][1] = cal_sa(pop, offspring, region_to_sort_ind[k], m_region_sort[k].first, pops);	//calculate sa
			data_dds_sa[k][1] = cal_cdi(pop, offspring, region_to_sort_ind[k], m_region_sort[k].first, pops,pro);	//calculate cdi
			data_dds_sa1[k] = &data_dds_sa[k];
		}

		/*std::vector<std::vector<Real>> test(5, std::vector<Real>(2));
		std::vector<std::vector<Real>*> test_ptr(5);
		std::vector<int> test_rank(5);
		test = { {2,1},{1,0},{4,2},{3,3},{0,4} };
		for (int i = 0; i < 5; ++i) {
			test_ptr[i] = &test[i];
		}
		std::vector<OptimizeMode> test_mode = { OptimizeMode::Maximize , OptimizeMode::Maximize };
		nd_sort::filter_sort(test_ptr, test_rank, test_mode);*/

		//std::vector<int> rank_dds_sa1(size_region);
		std::vector<OptimizeMode> opt_mode({ OptimizeMode::Maximize, OptimizeMode::Maximize }); //max dds and max sa
		nd_sort::filter_sort(data_dds_sa1, rank_dds_sa, opt_mode);
		if (pops <= pops_F1) {
			/*record dds and sa*/
			for (int k = 0; k < size_region; ++k) {
				m_region_sort1[k].second = data_dds_sa[k][0];
				m_region_sort2[k].second = data_dds_sa[k][1];
			}
			/*record rank after non-dominance sort*/
			for (int k = 0; k < size_region; ++k) {
				m_region_sort3[k].second = rank_dds_sa[k];
			}
		}

	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::sort_region_predict(const std::vector<std::vector<Real>>& prediction, const Population& pop, std::vector<int>& rank_dds_sa, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind, int pops, int pops_F1,Problem *pro) {
		cal_dds(offspring, region_to_sort_ind,pro);	//calculate dds
		int size_region = m_region_sort.size();
		if (rank_dds_sa.size() != size_region)
			rank_dds_sa.resize(size_region);
		//observe cdi
		if (m_observe) {
			std::vector<std::vector<Real>> temp(3, std::vector<Real>(CAST_CONOP(pro)->numberObjectives()));
			m_observe_cdi.resize(size_region, temp);
		}

		//non-dominance sort of dds and sa
		std::vector<std::vector<Real>> data_dds_sa(size_region, std::vector<Real>(2));
		std::vector<std::vector<Real>*> data_dds_sa1(size_region);
		for (int k = 0; k < size_region; ++k) {
			data_dds_sa[k][0] = m_region_sort[k].second;
			//data_dds_sa[k][1] = cal_sa(pop, offspring, region_to_sort_ind[k], m_region_sort[k].first, pops);	//calculate sa
			data_dds_sa[k][1] = cal_cdi_predict(prediction, pop, offspring, region_to_sort_ind[k], m_region_sort[k].first, pops, m_observe_cdi[k],pro);	//calculate cdi
			data_dds_sa1[k] = &data_dds_sa[k];
		}

		//oberve dds and cdi
		if (m_observe) {
			/*cout.setf(ios::fixed);
			for (int i = 0; i < size_region; ++i) {
				std::cout << i << ": " << m_region_sort[i].first << "\tdd: " << setprecision(5) << data_dds_sa[i][0] << "\tcdi: " << setprecision(5) << data_dds_sa[i][1] << std::endl;
				auto box = kdtree->get_box(m_region_sort[i].first);
				for (int j = 0; j < box.size(); ++j) {
					std::cout << "(" << box[j].first << "\t" << box[j].second << ")" << std::endl;
				}
			}
			std::cout << "min_max:" << std::endl;
			for (int j = 0; j < m_min_max.size(); ++j) {
				std::cout << "(" << m_min_max[j].first << "\t" << m_min_max[j].second << ")" << std::endl;
			}
			std::cout << std::endl;*/
			m_observe = false;
		}

		std::vector<OptimizeMode> opt_mode({ OptimizeMode::kMaximize, OptimizeMode::kMaximize }); //max dds and max sa
		nd_sort::filterSort(data_dds_sa1, rank_dds_sa, opt_mode);
		if (pops <= pops_F1) {
			/*record dds and sa*/
			for (int k = 0; k < size_region; ++k) {
				m_region_sort1[k].second = data_dds_sa[k][0];
				m_region_sort2[k].second = data_dds_sa[k][1];
			}
			/*record rank after non-dominance sort*/
			for (int k = 0; k < size_region; ++k) {
				m_region_sort3[k].second = rank_dds_sa[k];
			}
		}

	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::sort_region_new(const std::vector<std::vector<Real>>& ref_point, std::vector<int>& rank_norm, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind, Problem *pro) {
		update_min_max_predict(ref_point,pro);
		int M = CAST_CONOP(pro)->numberObjectives();
		int size_region = m_region_sort.size();
		if (m_observe) {
			std::vector<std::vector<Real>> temp(2, std::vector<Real>(M));
			m_observe_norm.resize(size_region, temp);
		}
		Real min_norm = 1e10;
		for (int k = 0; k < size_region; ++k) {
			std::vector<Real> x_k(M);	// x_k represent the subspace
			for (int j = 0; j < M; ++j) {
				Real temp_obj = 0;
				for (int i = 0; i < region_to_sort_ind[k].size(); ++i) {
					temp_obj += offspring[region_to_sort_ind[k][i]].objective(j);
				}
				x_k[j] = temp_obj / region_to_sort_ind[k].size();
			}
			std::vector<Real> z_star(M);
			if (!find_reference_point(ref_point, x_k, z_star,pro)) {
				for (int j = 0; j < M; ++j) {
					z_star[j] = m_min_max[j].first;
				}
			}
			m_region_sort[k].second = euclideanDistance(x_k.begin(), x_k.end(), z_star.begin());	//record
			if (m_region_sort[k].second < min_norm)
				min_norm = m_region_sort[k].second;
			if (m_observe) {
				m_observe_norm[k][0] = z_star; //reference point
				m_observe_norm[k][1] = x_k;	  //subspace point
			}
		}
		if (rank_norm.size() != size_region)
			rank_norm.resize(size_region, -1);
		for (int k = 0; k < size_region; ++k) {
			if (fabs(m_region_sort[k].second - min_norm) < 1e-10)
				rank_norm[k] = 0;
		}

		if (m_observe) {
			for (int k = 0; k < size_region; ++k) {
				m_region_sort_norm[k].second = m_region_sort[k].second; //record norm between z_star and x_k
			}
			m_observe = false;
		}

	}

	//non-dominance sort of norm and angle (or distance)
	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::sort_region_norm_angle(const std::vector<std::vector<Real>>& ref_point, const Population& pop, std::vector<int>& rank_norm_angle, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind, int pops,Problem *pro) {
		update_min_max_predict(ref_point,pro);
		int size_region = m_region_sort.size();
		int M = CAST_CONOP(pro)->numberObjectives();
		if (rank_norm_angle.size() != size_region)
			rank_norm_angle.resize(size_region);
		//observe norm
		if (m_observe) {
			std::vector<std::vector<Real>> temp_norm(2, std::vector<Real>(M));
			std::vector<std::vector<Real>> temp_angle(3, std::vector<Real>(M));
			m_observe_norm.resize(size_region, temp_norm);
			m_observe_angle.resize(size_region, temp_angle);
		}
		//non-dominance sort of norm and angle (or distance)
		std::vector<std::vector<Real>> data_norm_angle(size_region, std::vector<Real>(2));
		std::vector<std::vector<Real>*> data_norm_angle1(size_region);
		for (int k = 0; k < size_region; ++k) {
			/*****calculate norm*****/
			std::vector<Real> x_k(M);	// x_k represent the subspace
			std::vector<Real> z_ref(M); // z_ref represent the reference point of the subspace
			for (int j = 0; j < M; ++j) {
				Real temp_obj = 0;
				for (int i = 0; i < region_to_sort_ind[k].size(); ++i) {
					temp_obj += offspring[region_to_sort_ind[k][i]].objective(j);
				}
				x_k[j] = temp_obj / region_to_sort_ind[k].size();
			}
			if (!find_reference_point(ref_point, x_k, z_ref,pro)) {
				for (int j = 0; j < M; ++j) {
					z_ref[j] = m_min_max[j].first;
				}
			}
			//record norm
			if (m_observe) {
				m_observe_norm[k][0] = z_ref; //reference point
				m_observe_norm[k][1] = x_k;	  //subspace point
				m_observe_angle[k][0] = x_k;  //subspace point
			}
			//normalization x_k and z_ref
			for (int j = 0; j < M; ++j) {
				x_k[j] = (x_k[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
				z_ref[j] = (z_ref[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
			}
			//calculate norm
			data_norm_angle[k][0] = euclideanDistance(x_k.begin(), x_k.end(), z_ref.begin());

			/*****calculate angle (or distance)*****/
			std::pair<int, int> angle_idx(-1, -1);
			std::vector<std::pair<Real, int>> distance(pops, std::pair<Real, int>(-1, -1));
			for (int i = 0; i < pops; ++i) {
				std::vector<Real> ind_obj = pop[i].objective();
				for (int j = 0; j < M; ++j)	//normalizate ind_obj, i.e., pop[i].objective()
					ind_obj[j] = (ind_obj[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
				if (Dominance::kDominant == objectiveCompare(x_k, ind_obj, CAST_CONOP(pro)->optimizeMode()))
					distance[i].first = 1e10;
				else {
					distance[i].first = euclideanDistance(x_k.begin(), x_k.end(), ind_obj.begin());
				}
				distance[i].second = i;
			}
			std::sort(distance.begin(), distance.end(), [](std::pair<Real, int>a, std::pair<Real, int>b) {return a.first < b.first; });
			angle_idx.first = distance[0].second;
			angle_idx.second = distance[1].second;
			std::vector<Real> x_first = pop[angle_idx.first].objective();
			std::vector<Real> x_second = pop[angle_idx.second].objective();
			//record angle
			if (m_observe) {
				m_observe_angle[k][1] = x_first;//the first nearest Solution
				m_observe_angle[k][2] = x_second;//the second nearest Solution
			}
			//normalization x_first and x_second
			for (int j = 0; j < M; ++j) {
				x_first[j] = (x_first[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
				x_second[j] = (x_second[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
			}
			//calculate angle
			Real point_product, first_norm, second_norm;//point product (x_k-x_first) and (x_p-x_second), modules of (x_k-x_first), modules of (x_p-x_second)
			point_product = 0;
			for (int j = 0; j < M; ++j) {
				point_product += (x_first[j] - x_k[j]) * (x_second[j] - x_k[j]);
			}
			first_norm = euclideanDistance(x_k.begin(), x_k.end(), x_first.begin());
			second_norm = euclideanDistance(x_k.begin(), x_k.end(), x_second.begin());
			Real cos_theta = point_product / (first_norm * second_norm);
			if ((fabs(cos_theta) - 1.0) > 1e-6)
				std::cout << "calculate angle is error!" << std::endl;
			data_norm_angle[k][1] = acos(cos_theta);

			data_norm_angle1[k] = &data_norm_angle[k];	//for nonDominated sorting

		}

		if (m_observe) {
			std::cout << "min-max: " << std::endl;
			for (int j = 0; j < m_min_max.size(); ++j)
				std::cout << std::setw(10) << m_min_max[j].first << std::setw(10) << m_min_max[j].second << std::endl;
			m_observe = false;
		}

		std::vector<OptimizeMode> opt_mode({ OptimizeMode::kMinimize, OptimizeMode::kMaximize }); //min norm and max angle
		nd_sort::filterSort(data_norm_angle1, rank_norm_angle, opt_mode);
		if (pops <= m_count_F1) {
			/*record norm and angle*/
			for (int k = 0; k < size_region; ++k) {
				m_region_sort1[k].second = data_norm_angle[k][0];	//norm
				m_region_sort2[k].second = data_norm_angle[k][1];	//angle
			}
			/*record rank after non-dominance sort*/
			for (int k = 0; k < size_region; ++k) {
				m_region_sort3[k].second = rank_norm_angle[k];
			}
			//test map value
			/*std::cout << "norm\tangle" << std::endl;
			Real rank_max1 = 0, rank_min1 = 1e10, rank_max2 = 0, rank_min2 = 1e10;
			for (int k = 0; k < size_region; ++k) {
				if (rank_max1 < m_region_sort1[k].second)
					rank_max1 = m_region_sort1[k].second;
				if (m_region_sort1[k].second < rank_min1)
					rank_min1 = m_region_sort1[k].second;
				if (rank_max2 < m_region_sort2[k].second)
					rank_max2 = m_region_sort2[k].second;
				if (m_region_sort2[k].second < rank_min2)
					rank_min2 = m_region_sort2[k].second;
			}
			for (int k = 0; k < size_region; ++k) {
				Real map_norm, map_angle;
				map_norm = (m_region_sort1[k].second - rank_min1) * 255.0 / (rank_max1 - rank_min1);
				map_angle = 255.0 - (m_region_sort2[k].second - rank_min2) * 255.0 / (rank_max2 - rank_min2);
				std::cout << std::setw(10) << map_norm << std::setw(10) << map_angle << std::endl;
			}*/
		}
	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::sort_dd(std::vector<Solution>& offspring) {
		m_D.resize(offspring.size(), std::vector<Real>(offspring.size()));

		int size = offspring.size();
		//initialize rank
		for (int i = 0; i < size; ++i) {
			offspring[i].setFitness(-1);
		}
		//caculate dominance matrix m_D
		for (int i = 0; i < size; ++i) {
			for (int j = i + 1; j < offspring.size(); ++j) {
				m_D[i][j] = dd(offspring[i], offspring[j]);
				m_D[j][i] = -m_D[i][j];
			}
		}
		//calculate rank of Solution
		std::vector<double> ratio_dd(offspring.size());
		int rank = -1;
		int count = 0;
		while (count < size) {
			rank++;
			double max_dd = 0;
			for (int i = 0; i < size; ++i) {
				if (offspring[i].fitness() != -1)
					continue;
				double val_Dominant = 0, val_Dominated = 0;
				for (int j = 0; j < size; ++j) {
					if (offspring[j].fitness() != -1)
						continue;
					if (m_D[i][j] > 0)
						val_Dominant += m_D[i][j];
					else
						val_Dominated += std::fabs(m_D[i][j]);
				}
				if (std::fabs(val_Dominated) < 1e-6) {
					ratio_dd[i] = 1e6;
					max_dd = 1e6;
				}
				else if (std::fabs(val_Dominant) < 1e-6)
					ratio_dd[i] = -1e6;
				else {
					ratio_dd[i] = val_Dominant / val_Dominated;
					if (max_dd < ratio_dd[i])
						max_dd = ratio_dd[i];
				}
			}
			for (int i = 0; i < size; ++i) {
				if (std::fabs(ratio_dd[i] - max_dd) < 1e-6) {
					offspring[i].setFitness(rank);
					count++;
					ratio_dd[i] = 0;
					for (int j = 0; j < offspring.size(); ++j) {
						m_D[i][j] = 0;
						m_D[j][i] = 0;
					}
				}
			}
		}
		//test
		std::cout << "sort_dd rank_max: " << rank << std::endl;
		/*for (int i = 0; i < size; ++i) {
			std::cout << "[" << i << "]\t" << offspring[i].fitness() << std::endl;
		}*/
		//observation my dominance method
		/*std::vector<std::vector<int>> v1(rank + 1);
		for (int i = 0; i <= rank; ++i) {
			v1[i].reserve(offspring.size());
		}
		std::cout << "my dominace method" << std::endl;
		for (int rank_ = 0; rank_ <= rank; ++rank_) {
			std::cout << "[" << rank_ << "]: ";
			for (int i = 0; i < offspring.size(); ++i) {
				if (offspring[i].fitness() == rank_) {
					v1[rank_].push_back(i);
					std::cout << i << " ";
				}
			}
			std::cout << std::endl;
		}*/
		//compare the first rank Solution
		/*for (int i = 0; i < v1[0].size(); ++i) {
			if (offspring[v1[0][i]].type() == -1) {
				offspring[v1[0][i]].set_type(1);
				m_dd_pop++;
			}
			else if (offspring[v1[0][i]].type() == 0) {
				offspring[v1[0][i]].set_type(2);
				m_both_pop++;
				m_nsga2_pop--;
			}
		}*/
	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::sort_dd_unify(std::vector<Solution>& offspring,Problem *pro) {
		m_D.resize(offspring.size(), std::vector<Real>(offspring.size()));
		int size = offspring.size();
		//initialize rank
		for (int i = 0; i < size; ++i) {
			offspring[i].setFitness(-1);
		}
		//caculate dominance matrix m_D
		for (int i = 0; i < size; ++i) {
			for (int j = i + 1; j < offspring.size(); ++j) {
				m_D[i][j] = dd(offspring[i], offspring[j],pro);
				m_D[j][i] = -m_D[i][j];
			}
		}
		//calculate rank of Solution
		std::vector<double> ratio_dd(offspring.size());
		for (int i = 0; i < size; ++i) {
			double val_Dominant = 0, val_Dominated = 0;
			for (int j = 0; j < size; ++j) {
				if (m_D[i][j] > 0)
					val_Dominant += m_D[i][j];
				else
					val_Dominated += std::fabs(m_D[i][j]);
			}
			if (std::fabs(val_Dominated) < 1e-6) {
				ratio_dd[i] = 1e6;
				//max_dd = 1e6;
			}
			else if (std::fabs(val_Dominant) < 1e-6)
				ratio_dd[i] = -1e6;
			else {
				ratio_dd[i] = val_Dominant / val_Dominated;
				/*if (max_dd < ratio_dd[i])
					max_dd = ratio_dd[i];*/
			}
		}
		int rank = -1;
		int count = 0;
		while (count < size) {
			rank++;
			double max_dd = -1e6;
			for (int i = 0; i < size; ++i) {
				if (offspring[i].fitness() == -1 && max_dd < ratio_dd[i])
					max_dd = ratio_dd[i];
			}
			for (int i = 0; i < size; ++i) {
				if (std::fabs(ratio_dd[i] - max_dd) < 1e-6) {
					offspring[i].setFitness(rank);
					count++;
					//ratio_dd[i] = 0;
				}
			}
		}
		//test
		//std::cout << "sort_dd_unify rank_max: " << rank << std::endl;
		/*for (int i = 0; i < size; ++i) {
			std::cout << "[" << i << "]\t" << ratio_dd[i] << "\t" << offspring[i].fitness() << std::endl;
		}*/
	}


	template<typename Population, typename Solution>
	Real KDT_MOEA<Population, Solution>::dd(Solution& ind1, Solution& ind2,Problem *pro) {
		int better, worse;
		better = worse = 0;
		int M = CAST_CONOP(pro)->numberObjectives();
		for (int j = 0; j < M; ++j) {
			if (ind1.objective(j) < ind2.objective(j))
				better++;
			else if (ind2.objective(j) < ind1.objective(j))
				worse++;
		}
		return (Real)(better - worse) / M;
	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::sort_1_k_dominance(std::vector<Solution>& offspring, int m_k, Problem *pro) {
		std::vector<std::vector<Real>*> objs;
		for (auto& i : offspring)
			objs.emplace_back(&i.objective());
		std::vector<int> rank;
		dominance_1_k(objs, rank, CAST_CONOP(pro)->optimizeMode(), m_k);
		for (size_t i = 0; i < offspring.size(); ++i)
			offspring[i].setFitness(rank[i]);
	}
	template<typename Population, typename Solution>
	int KDT_MOEA<Population, Solution>::dominance_1_k(const std::vector<std::vector<Real>*>& data, std::vector<int>& rank, const std::vector<OptimizeMode>& opt_mode, int m_k) {
		const std::size_t popsize(data.size());
		if (popsize == 0) return 0;
		if (rank.size() != popsize)
			rank.resize(popsize);
		std::vector<int> rank_(popsize);
		std::vector<int> count(popsize);
		std::vector<std::vector<int>> cset(popsize, std::vector<int>(popsize));

		for (int i = 0; i < popsize; i++)
			rank[i] = -1;

		for (int k = 0; k < popsize; k++) {
			for (int j = 0; j < popsize; j++) {
				if (k != j) {
					auto compare_result = objective_compare_1_k(*data[j], *data[k], opt_mode, m_k);
					if (compare_result == Dominance::kDominant) {//*data[j]>*data[k]
						rank_[k]++;
					}
					else if (compare_result == Dominance::kDominated) {//*data[k]>*data[j]
						cset[k][count[k]] = j;
						count[k]++;
					}
				}
			}
		}
		int m_curRank = 0;
		std::vector<int> rank2(popsize);
		std::vector<int> number;
		int number_;
		int number_sum = 0;
		while (1) {
			int stop_count = 0;
			int min_rank = 1e6;
			number_ = 0;
			for (int k = 0; k < popsize; k++) {
				rank2[k] = rank_[k];
				if (rank[k] == -1 && rank_[k] < min_rank)
					min_rank = rank_[k];
			}
			auto i = cset.begin();
			for (int k = 0; k < popsize; k++, ++i) {
				if (rank[k] == -1 && rank_[k] == min_rank) {
					rank[k] = m_curRank;
					number_++;
					for (int j = 0; j < count[k]; j++) {
						int id = (*i)[j];
						rank2[id]--;
						stop_count++;
					}
				}
			}
			if (number_) {
				number.push_back(number_);
				number_sum += number_;
			}
			for (int k = 0; k < popsize; k++)
				rank_[k] = rank2[k];
			m_curRank++;
			if (stop_count == 0) {
				if (number_sum != popsize)
					std::cout << "error" << std::endl;
				return m_curRank;
			}
		}
	}
	template<typename Population, typename Solution>
	Dominance KDT_MOEA<Population, Solution>::objective_compare_1_k(const std::vector<Real>& a, const std::vector<Real>& b, const std::vector<OptimizeMode>& mode, int m_k) {
		if (a.size() != b.size())
			return Dominance::kNonComparable;

		int better = 0, worse = 0, equal = 0;
		for (decltype(a.size()) i = 0; i < a.size(); ++i) {
			if (mode[i] == OptimizeMode::kMinimize) {
				if (a[i] < b[i])
					better++;
				else if (a[i] > b[i])
					worse++;
				else
					equal++;
				/*if (a[i] < b[i]) {
					if (worse > 0)
						return Dominance::Non_Dominated;
					else
						++better;
				}
				else if (a[i] > b[i]) {
					if (better > 0)
						return Dominance::Non_Dominated;
					else
						++worse;
				}
				else {
					++equal;
				}*/
			}
			else {
				if (a[i] > b[i])
					better++;
				else if (a[i] < b[i])
					worse++;
				else
					equal++;
				/*if (a[i] > b[i]) {
					if (worse > 0)
						return Dominance::Non_Dominated;
					else
						++better;
				}
				else if (a[i] < b[i]) {
					if (better > 0)
						return Dominance::Non_Dominated;
					else
						++worse;
				}
				else {
					++equal;
				}*/
			}
		}
		if (equal == a.size())
			return Dominance::kEqual;
		else if (better >= ((Real)(better + worse) / (m_k + 1)))
			return Dominance::kDominant;
		else if (worse >= ((Real)(worse + better) / (m_k + 1)))
			return Dominance::kDominated;
		else
			return Dominance::kNonDominated;
	}



	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::cal_dds(const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind,Problem *pro) {
		//compute domination value
		int size = m_region_sort.size();
		int M = CAST_CONOP(pro)->numberObjectives();
		std::vector<std::vector<Real>> D(size, std::vector<Real>(size));
		std::vector<std::pair<Real, Real>> sub1(M), sub2(M);
		std::vector<Real> domination_value(size, 0);
		Real domination_max = 0;
		for (int i = 0; i < size; ++i) {
			Real Dominant_val = 0, Dominated_val = 0;
			for (int j = 0; j < size; ++j) {
				D[i][i] = 0;
				if (i < j) {
					//D[i][j] = cal_D_ij((std::vector<std::pair<Real, Real>>)kdtree->get_box(m_region_sort[i].first), (std::vector<std::pair<Real, Real>>)kdtree->get_box(m_region_sort[j].first));
					get_sub(sub1, offspring, region_to_sort_ind, i);
					get_sub(sub2, offspring, region_to_sort_ind, j);
					D[i][j] = cal_D_ij(sub1, sub2,pro);
					D[j][i] = -D[i][j];
				}
				if (D[i][j] > 0)
					Dominant_val += D[i][j];
				else if (D[i][j] < 0)
					Dominated_val += fabs(D[i][j]);
			}
			if (Dominated_val < 1.0e-6) {	//the best region
				m_region_sort[i].second = 0;
				domination_value[i] = 1.0e6;
			}
			else if (Dominant_val < 1.0e-6) {	//the worst region
				m_region_sort[i].second = -1;
				domination_value[i] = -1.0e6;
			}
			else {
				domination_value[i] = (Real)Dominant_val / Dominated_val;
				m_region_sort[i].second = domination_value[i];
				if (domination_max < domination_value[i])
					domination_max = domination_value[i];
			}

		}
		for (int i = 0; i < size; ++i) {
			if (abs(m_region_sort[i].second) < 1.0e-6)
				m_region_sort[i].second = domination_max + 1;
		}
		////sort subspace according to domination value			 
		//for (int i = 0; i < size; ++i) {
		//	for (int j = i; j < size; ++j) {
		//		if (domination_value[i] < domination_value[j]) {
		//			Real temp1 = domination_value[i];
		//			domination_value[i] = domination_value[j];
		//			domination_value[j] = temp1;
		//			auto temp2 = m_region_sort[i];
		//			m_region_sort[i] = m_region_sort[j];
		//			m_region_sort[j] = temp2;
		//		}
		//	}
		//}
		//test
		/*std::cout << "domination: " << std::endl;
		for (int i = 0; i < size; ++i) {
		std::cout << m_region_sort[i].second << " ";
		}
		std::cout << std::endl;*/
	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::get_sub(std::vector<std::pair<Real, Real>>& sub, const std::vector<Solution>& offspring, const std::vector<std::vector<int>>& region_to_sort_ind, int k) {
		for (int j = 0; j < sub.size(); ++j) {
			sub[j].first = sub[j].second = offspring[region_to_sort_ind[k][0]].objective()[j];
			for (int i = 0; i < region_to_sort_ind[k].size(); ++i) {
				int s = region_to_sort_ind[k][i];
				if (sub[j].first > offspring[region_to_sort_ind[k][i]].objective()[j])
					sub[j].first = offspring[region_to_sort_ind[k][i]].objective()[j];
				if (sub[j].second < offspring[region_to_sort_ind[k][i]].objective()[j])
					sub[j].second = offspring[region_to_sort_ind[k][i]].objective()[j];
			}
		}
	}

	template<typename Population, typename Solution>
	Real KDT_MOEA<Population, Solution>::cal_D_ij(const std::vector<std::pair<Real, Real>>& sub1, const std::vector<std::pair<Real, Real>>& sub2,Problem *pro) {
		int M = CAST_CONOP(pro)->numberObjectives();
		Real radio;
		Real better = 0, worse = 0;
		for (int i = 0; i < M; ++i) {
			if (sub1[i].second < sub2[i].first)
				better++;
			else if (sub2[i].second < sub1[i].first)
				worse++;
			else if (abs(sub1[i].second - sub2[i].first) < 1.0e-6) {
				if (sub1[i].first < sub2[i].second)
					better++;
				else if (sub2[i].first < sub1[i].second)
					worse++;
			}
			else if (abs(sub2[i].second - sub1[i].first) < 1.0e-6) {
				if (sub1[i].first < sub2[i].second)
					better++;
				else if (sub2[i].first < sub1[i].second)
					worse++;
			}

		}
		/*if (better != 0 && worse == 0)
		radio = (better / M);
		else if (better == 0 && worse != 0)
		radio = -(worse / M);
		else
		radio = 0;*/
		radio = (better - worse) / M;
		return radio;
	}

	/*template<typename Population, typename Solution>
	Real KDT_MOEA<Population, Solution>::cal_sa(const Population &pop, int k, int pops) {
		vector<Real> v_sa(pops);
		Real smallest_angle;
		int M = global::ms_global->m_problem->objective_size();
		vector<Real> x_k(M);		//it center point of kth subspace
		std::vector<std::pair<Real, Real>> box = kdtree->get_box(k);
		for (int j = 0; j < M; ++j) {
			x_k[j] = box[j].first + (box[j].second - box[j].first) / 2;
		}
		for (int i = 0; i < pops; ++i) {
			Real x_k_p, k_, p_;		//point product x_k and x_p, modules of x_k, modules of x_p
			x_k_p = 0;
			vector<Real> x_p = pop[i].objective();
			for (int j = 0; j < M; ++j) {
				x_k_p += x_k[j] * x_p[j];
			}
			vector<Real> x_0(M);
			k_ = euclidean_distance(x_k.begin(), x_k.end(), x_0.begin());
			p_ = euclidean_distance(x_p.begin(), x_p.end(), x_0.begin());
			v_sa[i] = x_k_p / (k_*p_);
			v_sa[i] = acos(v_sa[i]);
		}
		smallest_angle = v_sa[0];
		for (int i = 1; i < pops; ++i) {
			if (v_sa[i] < smallest_angle)
				smallest_angle = v_sa[i];
		}
		return smallest_angle;	//return smallest angle
	}*/

	template<typename Population, typename Solution>
	Real KDT_MOEA<Population, Solution>::cal_sa(const Population& pop, const std::vector<Solution>& offspring, const std::vector<int>& region_to_sort_ind_k, int k, int pops,Problem *pro) {
		std::vector<Real> v_sa(pops);
		Real smallest_angle;
		int M = pro->numberObjectives();
		std::vector<Real> x_k(M);	// x_k represent the subspace		
		std::vector<std::pair<Real, Real>> box = kdtree->getBox(k);
		/*for (int j = 0; j < M; ++j) {	//use center point of the subspace to represent this subspace
			x_k[j] = box[j].first + (box[j].second - box[j].first) / 2;
		}*/
		for (int j = 0; j < M; ++j) {	//use Solutions of the subspace to represent this subspace
			Real temp_obj = 0;
			for (int i = 0; i < region_to_sort_ind_k.size(); ++i) {
				temp_obj += offspring[region_to_sort_ind_k[i]].objective(j);
			}
			x_k[j] = temp_obj / region_to_sort_ind_k.size();
		}
		for (int j = 0; j < M; ++j) {	//normalization
			x_k[j] = (x_k[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
		}
		std::vector<Real> z_star(M);	// the MOP is considered to be minimal optimization, z_star is z* approximately
		/*for (int j = 0; j < M; ++j) {
			z_star[j] = m_min_max[j].first;
		}*/
		for (int j = 0; j < M; ++j) {	//normalization
			z_star[j] = 0;
		}

		for (int i = 0; i < pops; ++i) {
			Real x_k_p, k_, p_;		//point product (x_k-z*) and (x_p-z*), modules of (x_k-z*), modules of (x_p-z*)
			x_k_p = 0;
			std::vector<Real> x_p = pop[i].objective();
			for (int j = 0; j < M; ++j) {	//normalization
				x_p[j] = (x_p[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
			}
			for (int j = 0; j < M; ++j) {
				x_k_p += (x_k[j] - z_star[j]) * (x_p[j] - z_star[j]);
			}
			k_ = euclideanDistance(x_k.begin(), x_k.end(), z_star.begin());
			p_ = euclideanDistance(x_p.begin(), x_p.end(), z_star.begin());
			v_sa[i] = x_k_p / (k_ * p_);
			v_sa[i] = acos(v_sa[i]);
		}
		smallest_angle = v_sa[0];
		for (int i = 1; i < pops; ++i) {
			if (v_sa[i] < smallest_angle)
				smallest_angle = v_sa[i];
		}
		return smallest_angle;	//return smallest angle
	}

	template<typename Population, typename Solution>
	Real KDT_MOEA<Population, Solution>::cal_cdi(const Population& pop, const std::vector<Solution>& offspring, const std::vector<int>& region_to_sort_ind_k, int k, int pops,Problem *pro) {
		Real cdi;
		Real smallest_angle;
		Real normalizated_norm = 0;
		std::vector<Real> v_sa(pops);
		int M = CAST_CONOP(pro)->numberObjectives();
		std::vector<Real> x_k(M);	// x_k represent the subspace		
		std::vector<std::pair<Real, Real>> box = kdtree->getBox(k);
		/*for (int j = 0; j < M; ++j) {	//use center point of the subspace to represent this subspace
			x_k[j] = box[j].first + (box[j].second - box[j].first) / 2;
		}*/
		for (int j = 0; j < M; ++j) {	//use Solutions of the subspace to represent this subspace
			Real temp_obj = 0;
			for (int i = 0; i < region_to_sort_ind_k.size(); ++i) {
				temp_obj += offspring[region_to_sort_ind_k[i]].objective(j);
			}
			x_k[j] = temp_obj / region_to_sort_ind_k.size();
		}
		for (int j = 0; j < M; ++j) {	//normalization
			x_k[j] = (x_k[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
		}
		for (int j = 0; j < M; ++j) {	//calculate the norm of subspace(x_k)
			normalizated_norm += std::pow(x_k[j], 2);
		}
		normalizated_norm = std::sqrt(normalizated_norm);
		std::vector<Real> z_star(M);	// the MOP is considered to be minimal optimization, z_star is z* approximately
		/*for (int j = 0; j < M; ++j) {
			z_star[j] = m_min_max[j].first;
		}*/
		for (int j = 0; j < M; ++j) {	//normalization
			z_star[j] = 0;
		}

		for (int i = 0; i < pops; ++i) {
			Real x_k_p, k_, p_;		//point product (x_k-z*) and (x_p-z*), modules of (x_k-z*), modules of (x_p-z*)
			x_k_p = 0;
			std::vector<Real> x_p = pop[i].objective();
			for (int j = 0; j < M; ++j) {	//normalization
				x_p[j] = (x_p[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
			}
			for (int j = 0; j < M; ++j) {
				x_k_p += (x_k[j] - z_star[j]) * (x_p[j] - z_star[j]);
			}
			k_ = euclideanDistance(x_k.begin(), x_k.end(), z_star.begin());
			p_ = euclideanDistance(x_p.begin(), x_p.end(), z_star.begin());
			v_sa[i] = x_k_p / (k_ * p_);
			v_sa[i] = acos(v_sa[i]);
		}
		smallest_angle = v_sa[0];
		for (int i = 1; i < pops; ++i) {
			if (v_sa[i] < smallest_angle)
				smallest_angle = v_sa[i];
		}
		cdi = smallest_angle / normalizated_norm;
		return cdi;	//return convergence diversity indicator (cdi)
	}

	template<typename Population, typename Solution>
	Real KDT_MOEA<Population, Solution>::cal_cdi_predict(const std::vector<std::vector<Real>>& prediction, const Population& pop, const std::vector<Solution>& offspring, const std::vector<int>& region_to_sort_ind_k, int k, int pops, std::vector<std::vector<Real>>& observe_cdi_k,Problem *pro) {
		update_min_max_predict(prediction,pro);
		Real cdi_prediction;
		Real smallest_angle;
		Real norm = 0;
		std::vector<Real> v_sa(pops);
		int M = CAST_CONOP(pro)->numberObjectives();
		std::vector<Real> x_k(M);	// x_k represent the subspace		
		std::vector<std::pair<Real, Real>> box = kdtree->getBox(k);
		/*for (int j = 0; j < M; ++j) {	//use center point of the subspace to represent this subspace
			x_k[j] = box[j].first + (box[j].second - box[j].first) / 2;
		}*/
		for (int j = 0; j < M; ++j) {	//use Solutions of the subspace to represent this subspace
			Real temp_obj = 0;
			for (int i = 0; i < region_to_sort_ind_k.size(); ++i) {
				temp_obj += offspring[region_to_sort_ind_k[i]].objective(j);
			}
			x_k[j] = temp_obj / region_to_sort_ind_k.size();
		}

		//find reference point z_star, if finding is failed, then let z_star=z^*
		std::vector<Real> z_star(M);
		if (!find_reference_point(prediction, x_k, z_star,pro)) {
			for (int j = 0; j < M; ++j) {
				z_star[j] = m_min_max[j].first;
			}
		}
		/*std::pair<Real, int> dis(1e10, -1); //here z_star is the nearest prediction point
		for (int i = 0; i < prediction.size(); ++i) {
			Real temp = euclidean_distance(x_k.begin(), x_k.end(), prediction[i].begin());
			if (temp < dis.first) {
				dis.first = temp;
				dis.second = i;
			}
		}
		z_star = prediction[dis.second];*/
		if (m_observe) {
			observe_cdi_k[0] = z_star;//reference point
			observe_cdi_k[1] = x_k;	//subspace point
		}

		for (int j = 0; j < M; ++j) {	//normalization
			x_k[j] = (x_k[j] - m_min_max_predict[j].first) / (m_min_max_predict[j].second - m_min_max_predict[j].first);
			z_star[j] = (z_star[j] - m_min_max_predict[j].first) / (m_min_max_predict[j].second - m_min_max_predict[j].first);
		}
		norm = euclideanDistance(x_k.begin(), x_k.end(), z_star.begin());

		for (int i = 0; i < pops; ++i) {
			Real x_k_p, k_, p_;		//point product (x_k-z*) and (x_p-z*), modules of (x_k-z*), modules of (x_p-z*)
			x_k_p = 0;
			std::vector<Real> x_p = pop[i].objective();
			//for (int j = 0; j < M; ++j) {	//normalization
			//	x_p[j] = (x_p[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
			//}
			for (int j = 0; j < M; ++j) {
				x_k_p += (x_k[j] - z_star[j]) * (x_p[j] - z_star[j]);
			}
			k_ = euclideanDistance(x_k.begin(), x_k.end(), z_star.begin());
			p_ = euclideanDistance(x_p.begin(), x_p.end(), z_star.begin());
			v_sa[i] = x_k_p / (k_ * p_);
			v_sa[i] = acos(v_sa[i]);
		}
		smallest_angle = v_sa[0];
		int index_sa = 0;
		for (int i = 1; i < pops; ++i) {
			if (v_sa[i] < smallest_angle) {
				smallest_angle = v_sa[i];
				index_sa = i;
			}
		}

		//observe cdi
		/*observe_cdi_k[0] = z_star;//reference point
		observe_cdi_k[1] = x_k;	//subspace pint */
		if (m_observe) {
			observe_cdi_k[2] = pop[index_sa].objective();
		}

		cdi_prediction = smallest_angle / norm;
		return cdi_prediction;	//return convergence diversity indicator (cdi)
	}

	template<typename Population, typename Solution>
	bool KDT_MOEA<Population, Solution>::find_reference_point(const std::vector<std::vector<Real>>& prediction, const std::vector<Real>& s_k, std::vector<Real>& z_star,Problem *pro) {
		int k = kdtree->getRegionIdx(s_k);
		//test whether the left locate in the subspace
		/*auto box = kdtree->get_box(k);
		std::vector<Real> left(box.size());
		for (int j = 0; j < box.size(); ++j) {
			left[j] = box[j].first;
		}
		int k_ = kdtree->get_regionIdx(left);*/

		//std::vector<bool> flag(prediction.size(), false);
		std::vector<int> inner_and_dominating;
		std::vector<int> inner;
		std::vector<int> dominating;
		for (int i = 0; i < prediction.size(); ++i) {
			bool r1 = (k == kdtree->getRegionIdx(prediction[i]));
			bool r2 = (Dominance::kDominant == objectiveCompare(prediction[i], s_k, CAST_CONOP(pro)->optimizeMode()));
			if (r1 && r2) {
				inner_and_dominating.push_back(i);
			}
			else if (r1) {
				inner.push_back(i);
			}
			else if (r2) {
				dominating.push_back(i);
			}
		}
		bool result = false;
		if (inner_and_dominating.size() == 1) {
			z_star = prediction[inner_and_dominating[0]];
			result = true;
		}
		else if (inner_and_dominating.size() > 1) {
			std::pair<Real, int> dis(1e10, -1);
			for (int i = 0; i < inner_and_dominating.size(); ++i) {
				Real temp = euclideanDistance(s_k.begin(), s_k.end(), prediction[inner_and_dominating[i]].begin());
				if (temp < dis.first) {
					dis.first = temp;
					dis.second = inner_and_dominating[i];
				}
			}
			z_star = prediction[dis.second];
			result = true;
		}
		else {
			if (dominating.size() != 0) {
				std::pair<Real, int> dis(1e10, -1);
				for (int i = 0; i < dominating.size(); ++i) {
					Real temp = 0;
					for (int j = 0; j < s_k.size(); ++j) {
						temp += std::pow((s_k[j] - prediction[dominating[i]][j]) / (m_min_max_predict[j].second - m_min_max_predict[j].first), 2);
					}
					temp = sqrt(temp);
					//Real temp = euclidean_distance(s_k.begin(), s_k.end(), prediction[dominating[i]].begin());
					if (temp < dis.first) {
						dis.first = temp;
						dis.second = dominating[i];
					}
				}
				z_star = prediction[dis.second];
				result = true;
			}
		}
		return result;
	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::update_min_max(const Population& pop,Problem *pro) {
		int M = CAST_CONOP(pro)->numberObjectives();
		if (m_min_max.size() != M)
			m_min_max.resize(M);
		int size = pop.size();
		for (int i = 0; i < M; ++i) {
			m_min_max[i].first = 1.0e14;
			for (int j = 0; j < size; ++j) {
				if (pop[j].objective()[i] < m_min_max[i].first) {
					m_min_max[i].first = pop[j].objective()[i];
				}
			}
		}
		for (int i = 0; i < M; ++i) {
			m_min_max[i].second = -1 * 1.0e14;
			for (int j = 0; j < size; ++j) {
				if (m_min_max[i].second < pop[j].objective()[i]) {
					m_min_max[i].second = pop[j].objective()[i];
				}
			}
		}
	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::update_min_max(const std::vector<Solution>& offspring,Problem *pro) {
		int M = CAST_CONOP(pro)->numberObjectives();
		if (m_min_max.size() != M)
			m_min_max.resize(M);
		int size = offspring.size();
		for (int i = 0; i < M; ++i) {
			m_min_max[i].first = 1.0e14;
			for (int j = 0; j < size; ++j) {
				if (offspring[j].objective()[i] < m_min_max[i].first) {
					m_min_max[i].first = offspring[j].objective()[i];
				}
			}
		}
		for (int i = 0; i < M; ++i) {
			m_min_max[i].second = -1 * 1.0e14;
			for (int j = 0; j < size; ++j) {
				if (m_min_max[i].second < offspring[j].objective()[i]) {
					m_min_max[i].second = offspring[j].objective()[i];
				}
			}
		}
	}

	template<typename Population, typename Solution>
	void KDT_MOEA<Population, Solution>::update_min_max_predict(const std::vector<std::vector<Real>>& prediction,Problem *pro) {
		m_min_max_predict = m_min_max;
		int M = CAST_CONOP(pro)->numberObjectives();
		for (int j = 0; j < M; ++j) {
			for (int i = 0; i < prediction.size(); ++i) {
				if (prediction[i][j] < m_min_max_predict[j].first)
					m_min_max_predict[j].first = prediction[i][j];
				if (m_min_max_predict[j].second < prediction[i][j])
					m_min_max_predict[j].second = prediction[i][j];
			}
		}
	}

}

#endif	// !OFEC_KDT_MOEA_H
