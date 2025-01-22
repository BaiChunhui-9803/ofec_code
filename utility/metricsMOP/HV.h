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
// Created: 11 Mar 2024


#ifndef OFEC_HV_H
#define OFEC_HV_H

#include "../../core/definition.h"
#include "../../utility/random/newran.h"
#include "../functional.h"
#include <vector>
#include <numeric>

namespace ofec {
	void update_S(std::pair<Real, std::vector<std::vector<Real>>>& v1, std::vector<std::pair<Real, std::vector<std::vector<Real>>>>& s) {
		int m = 0;
		if (!s.empty()) {
			int n = s.size();
			for (int i = 0; i < n; i++) {
				if (ifSame(v1.second, s[i].second)) {
					s[i].first += v1.first;
					m = 1;
					break;
				}
			}
		}
		if (m == 0)
			s.push_back(v1);
	}
	std::vector<std::vector<Real>> insertElement(std::vector<Real>& p, int k, std::vector<std::vector<Real>>& p1) {
		bool flag1 = 0;
		bool flag2 = 0;
		std::vector<std::vector<Real>> q1;
		if (!p1.empty()) {
			auto hp = p1[0];
			while (!p1.empty() && hp[k] < p[k]) {
				q1.push_back(hp);
				auto temp = p1;
				p1.clear();
				for (auto begin = temp.begin() + 1, end = temp.end(); begin != end; ++begin) {
					p1.push_back(*begin);
				}
				if (!p1.empty())
					hp = p1[0];
			}
		}
		q1.push_back(p);
		int m = p.size();
		while (!p1.empty()) {
			auto q = p1[0];
			for (int i = k; i < m; i++) {
				if (p[i] < q[i])
					flag1 = 1;
				else if (p[i] > q[i])
					flag2 = 1;
			}
			if (!(flag1 == 1 && flag2 == 0))
				q1.push_back(p1[0]);
			auto temp = p1;
			p1.clear();
			for (auto begin = temp.begin() + 1, end = temp.end(); begin != end; ++begin) {
				p1.push_back(*begin);
			}
		}
		return q1;
	}
	std::vector<std::pair<Real, std::vector<std::vector<Real>>>> slice(std::vector<std::vector<Real>>& obj, int j, std::vector<Real>ref_point) {
		std::vector<Real> p = obj[0];
		std::vector<std::vector<Real>> p1;
		for (auto begin = obj.begin() + 1, end = obj.end(); begin != end; ++begin) {
			p1.push_back(*begin);
		}
		std::vector<std::pair<Real, std::vector<std::vector<Real>>>> return_value;
		std::vector<std::vector<Real>> q1;
		while (!p1.empty()) {
			auto temp = insertElement(p, j + 1, q1);
			q1 = temp;
			auto p_ = p1[0];
			std::pair<Real, std::vector<std::vector<Real>>> ss;
			ss.first = fabs(p[j] - p_[j]);
			ss.second = q1;
			update_S(ss, return_value);
			p = p_;
			auto tmp = p1;
			p1.clear();
			for (auto begin = tmp.begin() + 1, end = tmp.end(); begin != end; ++begin) {
				p1.push_back(*begin);
			}
		}
		auto temp = insertElement(p, j + 1, q1);
		q1 = temp;
		std::pair<Real, std::vector<std::vector<Real>>> ss;
		ss.first = fabs(p[j] - ref_point[j]);
		ss.second = q1;
		update_S(ss, return_value);
		return return_value;
	}

	//all point are in [0,1]
	template<typename TPopulation>
	Real hyperVolume(const TPopulation& pop, const std::vector<Real>& ref_point, std::vector<OptimizeMode>& opt_mode , Random* rnd) {
		std::vector<std::vector<Real>> pop_temp;
		int m = pop.size();
		for (int i = 0; i < m; i++) {
			pop_temp.push_back(pop[i].objective());
		}
		Real score = 0.;
		int pop_size = pop_temp.size();
		int obj_size = pop_temp[0].size();
		std::vector<Real> fmin(obj_size, 0.);
		std::vector<Real> fmax(obj_size, 1.);
		////normalization
		//for (int i = 0; i < pop_size; i++) {
		//	for (int j = 0; j < obj_size; j++) {
		//		pop_temp[i][j] = (pop_temp[i][j] - fmin[0]) / (fmax[0] - fmin[0]) / 1.1;//1.1 is used to avoid nor get Real max and mini value
		//	}
		//}
		//set reference point
		//std::vector<T> ref_point(obj_size, 1.);
		if (pop_size == 0)
			return score;
		//exact HV value
		else if (obj_size < 4) {
			//sort according to the first objective value
			sortrows(pop_temp, 0);
			std::vector<std::pair<Real, std::vector<std::vector<Real>>>> S;
			S.push_back(std::make_pair(1., pop_temp));
			for (int i = 0; i < obj_size - 1; i++) {
				std::vector<std::pair<Real, std::vector<std::vector<Real>>>> S_;
				for (int j = 0; j < S.size(); j++) {
					auto s_temp = slice(S[j].second, i, ref_point);
					for (int k = 0; k < s_temp.size(); k++) {
						std::pair<Real, std::vector<std::vector<Real>>> ss;
						ss.first = s_temp[k].first * S[j].first;
						ss.second = s_temp[k].second;
						update_S(ss, S_);
					}
				}
				S = S_;
			}
			for (int i = 0; i < S.size(); i++) {
				auto temp = *S[i].second.begin();
				score = score + S[i].first * fabs(temp.back() - ref_point.back());
			}
			return score;
		}
		else {
			//if obj size is over 3
			size_t num_sample = 200000;
			return hypervolumeMontoCarlo(pop, ref_point, num_sample, opt_mode, rnd);
		}
	}

	template<typename T>
	Real hyperVolume(std::vector<std::vector<T>>& pop, std::vector<Real>& ref_point, std::vector<OptimizeMode> opt_mode, Random* rnd) {
		std::vector<std::vector<T>> pop_temp = pop;
		Real score = 0.;
		int pop_size = pop_temp.size();
		int obj_size = pop_temp[0].size();
		std::vector<Real> fmin(obj_size, 0.);
		std::vector<Real> fmax(obj_size, 1.);
		////normalization
		//for (int i = 0; i < pop_size; i++) {
		//	for (int j = 0; j < obj_size; j++) {
		//		pop_temp[i][j] = (pop_temp[i][j] - fmin[0]) / (fmax[0] - fmin[0]) / 1.1;//1.1 is used to avoid nor get Real max and mini value
		//	}
		//}
		//set reference point
		//std::vector<T> ref_point(obj_size, 1.);
		if (pop_size == 0)
			return score;
		//exact HV value
		else if (obj_size < 4) {
			//sort according to the first objective value
			sortrows(pop_temp, 0);
			std::vector<std::pair<Real, std::vector<std::vector<Real>>>> S;
			S.push_back(std::make_pair(1., pop_temp));
			for (int i = 0; i < obj_size - 1; i++) {
				std::vector<std::pair<Real, std::vector<std::vector<Real>>>> S_;
				for (int j = 0; j < S.size(); j++) {
					auto s_temp = slice(S[j].second, i, ref_point);
					for (int k = 0; k < s_temp.size(); k++) {
						std::pair<Real, std::vector<std::vector<Real>>> ss;
						ss.first = s_temp[k].first * S[j].first;
						ss.second = s_temp[k].second;
						update_S(ss, S_);
					}
				}
				S = S_;
			}
			for (int i = 0; i < S.size(); i++) {
				auto temp = *S[i].second.begin();
				score = score + S[i].first * fabs(temp.back() - ref_point.back());
			}
			return score;
		}
		else {
			//if obj size is over 3
			size_t num_sample = 200000;
			return hypervolumeMontoCarlo(pop_temp, ref_point, num_sample, opt_mode, rnd);
		}
	}

	/*calculate Hypervulume by MonteCarlo method*/
	template<typename TPopulation>
	Real hypervolumeMontoCarlo(const TPopulation& pop, const std::vector<Real>& ref_point, size_t nSample, std::vector<OptimizeMode>& opt_mode, Random* rnd) {
		Real HV = 0;
		//calculate lower point
		std::vector<Real> lower(ref_point.size());
		for (size_t dim = 0; dim < lower.size(); ++dim) {
			lower[dim] = pop[0].objective(dim);
			for (size_t i = 1; i < pop.size(); ++i) {
				if (pop[i].objective(dim) < lower[dim])
					lower[dim] = pop[i].objective(dim);
			}
		}
		//sample,and calculate number in HV		
		size_t nInside = 0;
		for (size_t i = 0; i < nSample; ++i) {
			std::vector<Real> point(ref_point.size());	//sample in [lower, ref_point]
			for (size_t i = 0; i < ref_point.size(); ++i) {
				point[i] = rnd->uniform.nextNonStd(lower[i], ref_point[i]);
			}
			bool isInside = false;						//judge whether the sampling point is in HV
			for (size_t i = 0; i < pop.size(); ++i) {
				if (objectiveCompare<>(pop[i].objective(), point, opt_mode) == Dominance::kDominant) {
					isInside = true;
					break;
				}
			}
			if (isInside)
				++nInside;
		}
		//calculate HV, HV/V = nInside/nSample , V=(ref1-low1)*(ref2-low2)*...
		Real V = std::inner_product(ref_point.begin(), ref_point.end(), lower.begin(), 1.0, std::multiplies<Real>(), std::minus<Real>());
		HV = nInside * V / nSample;
		return HV;
	}

	/*calculate Hypervulume by MonteCarlo method*/
	template<typename T>
	Real hypervolumeMontoCarlo(const std::vector<std::vector<T>>& pop, const std::vector<Real>& ref_point, size_t nSample, std::vector<OptimizeMode>& opt_mode, Random* rnd) {
		Real HV = 0;
		//calculate lower point
		std::vector<Real> lower(ref_point.size());
		for (size_t dim = 0; dim < lower.size(); ++dim) {
			lower[dim] = pop[0][dim];
			for (size_t i = 1; i < pop.size(); ++i) {
				if (pop[i][dim] < lower[dim])
					lower[dim] = pop[i][dim];
			}
		}
		//sample,and calculate number in HV		
		size_t nInside = 0;
		for (size_t i = 0; i < nSample; ++i) {
			std::vector<Real> point(ref_point.size());	//sample in [lower, ref_point]
			for (size_t i = 0; i < ref_point.size(); ++i) {
				point[i] = rnd->uniform.nextNonStd(lower[i], ref_point[i]);
			}
			bool isInside = false;						//judge whether the sampling point is in HV
			for (size_t i = 0; i < pop.size(); ++i) {
				if (objectiveCompare<>(pop[i], point, opt_mode) == Dominance::kDominant) {
					isInside = true;
					break;
				}
			}
			if (isInside)
				++nInside;
		}
		//calculate HV, HV/V = nInside/nSample , V=(ref1-low1)*(ref2-low2)*...
		Real V = std::inner_product(ref_point.begin(), ref_point.end(), lower.begin(), 1.0, std::multiplies<Real>(), std::minus<Real>());
		HV = nInside * V / nSample;
		return HV;
	}
}

#endif
