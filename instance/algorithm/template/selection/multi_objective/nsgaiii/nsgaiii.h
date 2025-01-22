/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia & Junchen Wang 
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 11 Jan 2015 by Changhe Li
// Last modified: 15 Aug 2019 by Junchen Wang (Email: wangjunchen.chris@gmail.com) 

/* --------------------------------NSGAIII-------------------------------------------      
 Deb and Jain, "An Evolutionary Many-Objective Optimization Algorithm Using
 Reference-point Based Non-dominated Sorting Approach, Part I: Solving Problems with
 Box Constraints," IEEE Transactions on Evolutionary Computation, to appear.
 http://dx.doi.org/10.1109/TEVC.2013.2281535
 ---------------------------------------------------------------------------------- */

#ifndef OFEC_NSGAIII_H
#define OFEC_NSGAIII_H

#include "../../../../../../utility/nondominated_sorting/fast_sort.h"
#include "reference_point.h"
#include <algorithm>

#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	template<typename TInd>
	class NSGAIII {
	public:
		NSGAIII(size_t size_pop, size_t size_obj, std::vector<OptimizeMode> opt_mode, Problem *pro) {
			setDefaultParam(pro);
			reference_point::generateRefPoints(&mv_rps, size_obj, mv_obj_division_p_);
			mvv_off_conv_obj.resize(2 * size_pop);
			for (int i = 0; i < 2 * size_pop; i++)
				mvv_off_conv_obj[i].resize(size_obj);
		}
		void survivorSelection(std::vector<std::unique_ptr<TInd>>& parent, std::vector<TInd>& offspring, Problem *pro, Random *rnd);
	private:
		void setDefaultParam(Problem *pro);
		void nondominatedSorting(std::vector<TInd>& offspring, Problem *pro);
		std::vector<Real> translateObjectives(const std::vector<std::vector<int> >& fronts, std::vector<TInd>& offspring, Problem *pro);
		void findExtremePoints(std::vector<size_t>* extreme_points, const std::vector<std::vector<int> >& fronts, std::vector<TInd>& offspring, Problem *pro);
		void constructHyperplane(std::vector<Real>* pintercepts, const std::vector<size_t>& extreme_points, std::vector<TInd>& offspring, Problem *pro);
		void normalizeObjectives(const std::vector<std::vector<int> >& fronts, const std::vector<Real>& intercepts, const std::vector<Real>& ideal_point, Problem *pro);
		size_t findNicheRefPoint(const std::vector<RefPoint>& rps, Random *rnd);
		int selectClusterMember(const RefPoint& rp, Random *rnd);
	protected:
		std::vector<size_t> mv_obj_division_p_;
		std::vector<RefPoint> mv_rps;
		std::vector<std::vector<Real>> mvv_off_conv_obj;
	};

	template<typename TInd>
	void NSGAIII<TInd>::setDefaultParam(Problem *pro) {
		//auto& v = pro;
		//if (v.name().find("DTLZ") != std::string::npos) {
		//	int num_obj = v.numberObjectives();
		//	if (num_obj == 3)
		//		mv_obj_division_p_.resize(1, 12);
		//	else if (num_obj == 5)
		//		mv_obj_division_p_.resize(1, 6);
		//	else if (num_obj == 8)
		//	{
		//		mv_obj_division_p_.resize(2);
		//		mv_obj_division_p_[0] = 3;
		//		mv_obj_division_p_[1] = 2;
		//	}
		//	return;
		//}
		mv_obj_division_p_.resize(1, 12);
	}

	template<typename TInd>
	void NSGAIII<TInd>::nondominatedSorting(std::vector<TInd>& offspring, Problem *pro) {
		std::vector<std::vector<Real>*> objs;
		for (auto& i : offspring)
			objs.emplace_back(&i.objective());
		std::vector<int> rank;
		nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(pro)->optimizeMode());
		for (size_t i = 0; i < offspring.size(); ++i)
			offspring[i].setFitness(rank[i]);
	}

	template<typename TInd>
	void NSGAIII<TInd>::survivorSelection(std::vector<std::unique_ptr<TInd>>& parent, std::vector<TInd>& offspring, Problem *pro, Random *rnd) {
		std::vector<RefPoint> rps = mv_rps;
		std::vector<std::vector<int> > fronts;
		int rank = 0;
		int count = 0;
		int size = offspring.size();
		// ---------- Step 4 in Algorithm 1: non-dominated sorting ----------
		nondominatedSorting(offspring, pro);
		while (1)
		{
			std::vector<int> temp;
			for (int i = 0; i < size; i++)
			{
				if (offspring[i].fitness() == rank)
				{
					temp.push_back(i);
					++count;
				}
			}
			fronts.push_back(temp);
			if (count == size) break;
			++rank;
		}

		// ---------- Steps 5-7 in Algorithm 1 ----------
		std::vector<size_t> considered; // St
		int last = 0, next_size = 0;
		while (next_size < parent.size())
		{
			next_size += fronts[last].size();
			last += 1;
		}
		fronts.erase(fronts.begin() + last, fronts.end()); // remove useless Solutions

		count = 0;
		for (size_t t = 0; t < fronts.size() - 1; t += 1)
			for (size_t i = 0; i < fronts[t].size(); i += 1)
				*(parent[count++]) = offspring[fronts[t][i]];

		// ---------- Steps 9-10 in Algorithm 1 ----------
		if (count == parent.size()) return;


		// ---------- Step 14 / Algorithm 2 ----------
		std::vector<Real> ideal_point = translateObjectives(fronts, offspring, pro);

		std::vector<size_t> extreme_points;
		findExtremePoints(&extreme_points, fronts, offspring, pro);

		std::vector<Real> intercepts;
		constructHyperplane(&intercepts, extreme_points, offspring, pro);

		normalizeObjectives(fronts, intercepts, ideal_point, pro);

		// ---------- Step 15 / Algorithm 3, Step 16 ----------
		reference_point::associate(&rps, mvv_off_conv_obj, fronts);

		// ---------- Step 17 / Algorithm 4 ----------
		while (count < parent.size())
		{
			size_t min_rp = findNicheRefPoint(rps, rnd);

			int chosen = selectClusterMember(rps[min_rp], rnd);
			if (chosen < 0) // no potential member in Fl, disregard this reference point
			{
				rps.erase(rps.begin() + min_rp);
			}
			else
			{
				rps[min_rp].addMember();
				rps[min_rp].removePotentialMember(chosen);
				*(parent[count++]) = offspring[chosen];
			}
		}
	}

	template<typename TInd>
	std::vector<Real> NSGAIII<TInd>::translateObjectives(const std::vector<std::vector<int>>& fronts, std::vector<TInd>& offspring, Problem *pro) {
		int numObj = pro->numberObjectives();
		std::vector<Real> ideal_point(numObj);

		for (int f = 0; f < numObj; f += 1)
		{
			Real minf = std::numeric_limits<Real>::max();
			for (size_t i = 0; i < fronts[0].size(); i += 1) // min values must appear in the first front
			{
				minf = std::min(minf, offspring[fronts[0][i]].objective()[f]);
			}
			ideal_point[f] = minf;

			for (size_t t = 0; t < fronts.size(); t += 1)
			{
				for (size_t i = 0; i < fronts[t].size(); i += 1)
				{
					size_t ind = fronts[t][i];
					mvv_off_conv_obj[ind][f] = offspring[ind].objective()[f] - minf;
				}
			}
		}

		return ideal_point;
	}

	template<typename TInd>
	void NSGAIII<TInd>::findExtremePoints(std::vector<size_t>* extreme_points, const std::vector<std::vector<int>>& fronts, std::vector<TInd>& offspring, Problem *pro) {
		int numObj = pro->numberObjectives();
		std::vector<size_t>& exp = *extreme_points;
		exp.clear();

		for (size_t f = 0; f < numObj; f += 1)
		{
			std::vector<Real> w(numObj, 0.000001);
			w[f] = 1.0;

			Real min_ASF = std::numeric_limits<Real>::max();
			size_t min_indv = fronts[0].size();

			for (size_t i = 0; i < fronts[0].size(); i += 1)  // only consider the Solutions in the first front
			{
				Real asf = math_aux::ASF(offspring[fronts[0][i]].objective(), w);
				if (asf < min_ASF)
				{
					min_ASF = asf;
					min_indv = fronts[0][i];
				}
			}

			exp.push_back(min_indv);
		}
	}

	template<typename TInd>
	void NSGAIII<TInd>::constructHyperplane(std::vector<Real>* pintercepts, const std::vector<size_t>& extreme_points, std::vector<TInd>& offspring, Problem *pro) {
		// Check whether there are duplicate extreme points.
	// This might happen but the original paper does not mention how to deal with it.
		int numObj = pro->numberObjectives();
		bool duplicate = false;
		for (size_t i = 0; !duplicate && i < extreme_points.size(); i += 1)
		{
			for (size_t j = i + 1; !duplicate && j < extreme_points.size(); j += 1)
			{
				duplicate = (extreme_points[i] == extreme_points[j]);
			}
		}

		std::vector<Real>& intercepts = *pintercepts;
		intercepts.assign(numObj, 0);

		if (duplicate) // cannot construct the unique hyperplane (this is a casual method to deal with the condition)
		{
			for (size_t f = 0; f < intercepts.size(); f += 1)
			{
				// extreme_points[f] stands for the Solution with the largest value of objective f
				intercepts[f] = offspring[extreme_points[f]].objective()[f];
			}
		}
		else
		{
			// Find the equation of the hyperplane
			std::vector<Real> b(numObj, 1.0);
			std::vector<std::vector<Real>> A;
			for (size_t p = 0; p < extreme_points.size(); p += 1)
			{
				A.push_back(offspring[extreme_points[p]].objective());
			}
			std::vector<Real> x;
			math_aux::guassianElimination(&x, A, b);

			// Find intercepts
			for (size_t f = 0; f < intercepts.size(); f += 1)
			{
				intercepts[f] = 1.0 / x[f];
			}
		}
	}

	template<typename TInd>
	void NSGAIII<TInd>::normalizeObjectives(const std::vector<std::vector<int>>& fronts, const std::vector<Real>& intercepts, const std::vector<Real>& ideal_point, Problem *pro) {
		int numObj = pro->numberObjectives();
		for (size_t t = 0; t < fronts.size(); t += 1)
		{
			for (size_t i = 0; i < fronts[t].size(); i += 1)
			{
				size_t ind = fronts[t][i];
				for (size_t f = 0; f < numObj; f += 1)
				{
					if (fabs(intercepts[f] - ideal_point[f]) > 10e-10) // avoid the divide-by-zero error
						mvv_off_conv_obj[ind][f] = mvv_off_conv_obj[ind][f] / (intercepts[f] - ideal_point[f]);
					else
						mvv_off_conv_obj[ind][f] = mvv_off_conv_obj[ind][f] / 10e-10;
				}
			}
		}
	}
	template<typename TInd>
	size_t NSGAIII<TInd>::findNicheRefPoint(const std::vector<RefPoint>& rps, Random *rnd) {
		// find the minimal cluster size
		size_t min_size = std::numeric_limits<size_t>::max();
		for (size_t r = 0; r < rps.size(); r += 1)
		{
			min_size = std::min(min_size, rps[r].memberSize());
		}

		// find the reference points with the minimal cluster size Jmin
		std::vector<size_t> min_rps;
		for (size_t r = 0; r < rps.size(); r += 1)
		{
			if (rps[r].memberSize() == min_size)
			{
				min_rps.push_back(r);
			}
		}

		// return a random reference point (j-bar)
		auto r = rnd->uniform.nextNonStd((size_t)0, min_rps.size());
		return min_rps[r];
	}

	template<typename TInd>
	int NSGAIII<TInd>::selectClusterMember(const RefPoint& rp, Random *rnd) {
		int chosen = -1;
		if (rp.hasPotentialMember())
		{
			if (rp.memberSize() == 0) // currently has no member
			{
				chosen = rp.findClosestMember();
			}
			else
			{
				chosen = rp.randomMember(rnd);
			}
		}
		return chosen;
	}
}



#endif // !OFEC_NSGAIII_H
