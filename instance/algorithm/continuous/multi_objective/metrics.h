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
// Created: 18 November 2020


#ifndef OFEC_METRICS_H
#define OFEC_METRICS_H

#include "../../../../core/definition.h"
#include "../../../../core/global.h"
#include "../../../../core/algorithm/population.h"
#include "../../../../core/problem/solution.h"
#include "../../../../utility/functional.h"
#include <vector>
#include <numeric>

namespace ofec {
	template<typename T>
	void updateS(std::pair<T, std::vector<std::vector<T>>>& v1, std::vector<std::pair<T, std::vector<std::vector<T>>>>& s) {
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

	template<typename T>
	std::vector<std::vector<T>> insertElement(std::vector<T>& p, int k, std::vector<std::vector<T>>& p1) {
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

	template<typename T>
	std::vector<std::pair<T, std::vector<std::vector<T>>>> slice(std::vector<std::vector<T>>& obj, int j, std::vector<T>ref_point) {
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
			updateS(ss, return_value);
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
		updateS(ss, return_value);
		return return_value;
	}

	/*calculate Hypervulume by MonteCarlo method*/
	template<typename Ind>
	Real hvMontoCarloPop(const Population<Ind>& pop, const std::vector<Real>& ref_point, size_t nSample, Problem *pro, Random *rnd) {
		Real HV = 0;
		//calculate lower point
		std::vector<Real> lower(ref_point.size());
		for (size_t dim = 0; dim < lower.size(); ++dim) {
			lower[dim] = pop[0].objective()[dim];
			for (size_t i = 1; i < pop.size(); ++i) {
				if (pop[i].objective()[dim] < lower[dim])
					lower[dim] = pop[i].objective()[dim];
			}
		}
		//sample,and calculate number in HV
		//std::vector<Real> lower(ref_point.size(), 0.);
		size_t nInside = 0;
		for (size_t i = 0; i < nSample; ++i) {
			std::vector<Real> point(ref_point.size());	//sample in [lower, ref_point]
			for (size_t i = 0; i < ref_point.size(); ++i) {
				point[i] = rnd->uniform.nextNonStd(lower[i], ref_point[i]);
			}
			bool isInside = false;						//judge whether the sampling point is in HV
			for (size_t i = 0; i < pop.size(); ++i) {
				if (objectiveCompare<>(pop[i].objective(), point, CAST_CONOP(pro)->optimizeMode()) == Dominance::kDominant) {
					isInside = true;
					break;
				}
			}
			if (isInside)
				++nInside;
		}
		//calculate HV, HV/V = nInside/nSample , V=(ref1-low1)*(ref2-low2)*...
		//Real V = std::inner_product(ref_point.begin(), ref_point.end(), lower.begin(), 1.0, std::multiplies<Real>(), std::minus<Real>());
		Real V = 1.;
		for (size_t i = 0; i < ref_point.size(); ++i) {
			V *= (ref_point[i] - lower[i]);
		}
		HV = nInside * V / nSample;
		return HV;
	}

	/*calculate Hypervulume by MonteCarlo method*/
	template<typename T>
	Real hvMontoCarloVector(const std::vector<std::vector<T>>& pop, const std::vector<T>& ref_point, size_t nSample, Problem *pro, Random *rnd) {
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
		//std::vector<Real> lower(ref_point.size(), 0.);
		size_t nInside = 0;
		for (size_t i = 0; i < nSample; ++i) {
			std::vector<Real> point(ref_point.size());	//sample in [lower, ref_point]
			for (size_t i = 0; i < ref_point.size(); ++i) {
				point[i] = rnd->uniform.nextNonStd(lower[i], ref_point[i]);
			}
			bool isInside = false;						//judge whether the sampling point is in HV
			for (size_t i = 0; i < pop.size(); ++i) {
				if (objectiveCompare<>(pop[i], point, CAST_CONOP(pro)->optimizeMode()) == Dominance::kDominant) {
					isInside = true;
					break;
				}
			}
			if (isInside)
				++nInside;
		}
		//calculate HV, HV/V = nInside/nSample , V=(ref1-low1)*(ref2-low2)*...
		//Real V = std::inner_product(ref_point.begin(), ref_point.end(), lower.begin(), 1.0, std::multiplies<Real>(), std::minus<Real>());
		Real V = 1.;
		for (size_t i = 0; i < ref_point.size(); ++i) {
			V *= (ref_point[i] - lower[i]);
		}
		HV = nInside * V / nSample;
		return HV;
	}

	//all point are in [0,1], ref_point is the maximize value of each obj.
	//ref_point is from PF samples
	template<typename Ind>
	Real hypervolumePop(const Population<Ind>& pop, const std::vector<std::vector<Real>>& ref_points,Problem *pro,Random *rnd) {
		std::vector<std::vector<Real>> pop_temp;
		int m = pop.size();
		for (int i = 0; i < m; i++) {
			pop_temp.push_back(pop[i].objective());
		}
		Real score = 0.;
		int pop_size = pop_temp.size();
		int obj_size = pop_temp[0].size();
		std::vector<Real> fmin(obj_size, INT16_MAX);
		std::vector<Real> fmax(obj_size, -1. * INT16_MAX);

		for (int i = 0; i < pop_size; i++) {
			for (int j = 0; j < obj_size; j++) {
				if (pop_temp[i][j] < fmin[j]) {
					fmin[j] = pop_temp[i][j];
				}
			}
		}
		for (int i = 0; i < fmin.size(); ++i) {
			if (fmin[i] > 0)
				fmin[i] = 0;
		}

		for (int i = 0; i < ref_points.size(); i++) {
			for (int j = 0; j < obj_size; j++) {
				if (ref_points[i][j] > fmax[j]) {
					fmax[j] = ref_points[i][j];
				}
			}
		}

		//normalization
		for (int i = 0; i < pop_size; i++) {
			for (int j = 0; j < obj_size; j++) {
				pop_temp[i][j] = (pop_temp[i][j] - fmin[j]) / (fmax[j] - fmin[j]) / 1.1;//1.1 is used to avoid nor get Real max and mini value
			}
		}
		//set reference point
		std::vector<Real> ref_point(obj_size, 1.);
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
						updateS(ss, S_);
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
			size_t num_sample = 1000000;
			return hvMontoCarloVector(pop_temp, ref_point, num_sample,pro,rnd);
		}
	}

    template<typename T>
	Real hypervolumeVector(std::vector<std::vector<T>>& pop, const std::vector<std::vector<T>>& ref_points,Problem *pro,Random *rnd) {
		std::vector<std::vector<T>> pop_temp = pop;
		Real score = 0.;
		int pop_size = pop_temp.size();
		int obj_size = pop_temp[0].size();
		std::vector<Real> fmin(obj_size, INT16_MAX);
		std::vector<Real> fmax(obj_size, -1. * INT16_MAX);
		for (int i = 0; i < pop_size; i++) {
			for (int j = 0; j < obj_size; j++) {
				if (pop[i][j] < fmin[j]) {
					fmin[j] = pop[i][j];
				}
			}
		}
		for (int i = 0; i < fmin.size(); ++i) {
			if (fmin[i] > 0)
				fmin[i] = 0;
		}

		for (int i = 0; i < ref_points.size(); i++) {
			for (int j = 0; j < obj_size; j++) {
				if (ref_points[i][j] > fmax[j]) {
					fmax[j] = ref_points[i][j];
				}
			}
		}

		//normalization
		for (int i = 0; i < pop_size; i++) {
			for (int j = 0; j < obj_size; j++) {
				pop_temp[i][j] = (pop_temp[i][j] - fmin[j]) / (fmax[j] - fmin[j]) / 1.1;//1.1 is used to avoid nor get Real max and mini value
			}
		}
		//delete items over one
		std::vector<std::vector<Real>> pop_obj;
		for (size_t i = 0; i < pop_temp.size(); ++i) {
			bool flag = false;
			for (size_t j = 0; j < pop_temp[i].size(); ++j) {
				if (pop_temp[i][j] > 1) {
					flag = true;
					break;
				}
			}
			if (!flag) {
				pop_obj.emplace_back(pop_temp[i]);
			}
		}
		//set reference point
		std::vector<T> ref_point(obj_size, 1.);
		if (pop_size == 0)
			return score;
		//exact HV value
		else if (obj_size < 4) {
			//sort according to the first objective value
			sortrows(pop_obj, 0);
			std::vector<std::pair<Real, std::vector<std::vector<Real>>>> S;
			S.emplace_back(std::make_pair(1., pop_obj));
			for (int i = 0; i < obj_size - 1; i++) {
				std::vector<std::pair<Real, std::vector<std::vector<Real>>>> S_;
				for (int j = 0; j < S.size(); j++) {
					auto s_temp = slice(S[j].second, i, ref_point);
					for (int k = 0; k < s_temp.size(); k++) {
						std::pair<Real, std::vector<std::vector<Real>>> ss;
						ss.first = s_temp[k].first * S[j].first;
						ss.second = s_temp[k].second;
						updateS(ss, S_);
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
			size_t num_sample = 1000000;
			return hvMontoCarloVector(pop_obj, ref_point, num_sample,pro,rnd);
		}
	}
}

#endif // !OFEC_METRICS_H
