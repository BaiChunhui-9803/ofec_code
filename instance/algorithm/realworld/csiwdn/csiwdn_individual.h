/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*********************************************************************************/
// updated Apr 4, 2018 by Li Zhou

#ifndef OFEC_CSIWDN_Solution_H
#define OFEC_CSIWDN_Solution_H

#include "../../../../core/problem/solution.h"
#include "../../../../core/environment/environment.h"
#include "../../../problem/realworld/csiwdn/csiwdn_encoding.h"

#include <cmath>
#include <fstream>
#include <cstring>

namespace ofec {

	class IndCSIWDN : public Solution<VarCSIWDN>	{
	protected:
		SolutionType m_pv, m_pu;
		bool m_improved = false;
	public:
		bool isImproved() const{
			return m_improved;
		}
		IndCSIWDN(size_t number_objectives,size_t num_cons,size_t dim) : 
			Solution(number_objectives, num_cons, dim),
			m_pv(number_objectives, num_cons, dim), 
			m_pu(number_objectives, num_cons, dim) {}

		IndCSIWDN(const IndCSIWDN& indi) : Solution(indi), m_pv(indi.m_pv), m_pu(indi.m_pu) {}
		IndCSIWDN(IndCSIWDN&& indi) : Solution(std::move(indi)), m_pv(std::move(indi.m_pv)), m_pu(std::move(indi.m_pu)) {}
		IndCSIWDN(Solution<VarCSIWDN>&& indi) :Solution(std::move(indi)), m_pv(indi), m_pu(indi) {};
		IndCSIWDN(const Solution<VarCSIWDN>& indi) :Solution(indi), m_pv(indi), m_pu(indi) {};

		IndCSIWDN& operator=(IndCSIWDN& indi);
		IndCSIWDN& operator=(IndCSIWDN&& indi);


		void mutateFirstPart(const std::vector<std::vector<Real>> &prob, Environment *env, Random *rnd, const std::pair<int, int>& source_index);
		void mutateSecondPart(Environment *env,Real F, int jdx, const std::vector<std::vector<Real>> &prob,
			Solution<VarCSIWDN> *r1,
			Solution<VarCSIWDN> *r2,
			Solution<VarCSIWDN> *r3,
			Solution<VarCSIWDN> *r4 = 0,
			Solution<VarCSIWDN> *r5 = 0);
		SolutionType& mpu() {
			return m_pu;
		}
		SolutionType& mpv() {
			return m_pv;
		}

		void recombine(Environment *env, Random *rnd, size_t jdx, Real CR);
		int select(Environment *env, bool is_stable, const std::pair<int, int>& source_index);
		void coverFirstPart(VarCSIWDN& indi, Environment *env, const std::pair<int,int> &source_index);
		void coverSecondPart(VarCSIWDN& indi, Environment *env, const std::pair<int, int>& source_index);
		bool sameLocation(Environment *env, IndCSIWDN & indi, const std::pair<int, int>& source_index);
		void initialize(Environment *env, Random *rnd);
		
		SolutionType& trial();
	protected:
		void mutateMass(Environment *env, Real F, int jdx, Solution<VarCSIWDN> *r1,
			Solution<VarCSIWDN> *r2,
			Solution<VarCSIWDN> *r3,
			Solution<VarCSIWDN> *r4 = 0,
			Solution<VarCSIWDN> *r5 = 0);
		void mutateLocation(Environment *env, Random *rnd, const std::vector<std::vector<Real>>& prob, const std::pair<int, int>& source_index);
		void mutateDuration(Environment *env, const std::vector<std::vector<Real>> &prob, const std::pair<int, int>& source_index);

	public:
		Real m_distance_fitness;
		Real m_pu_distance_fitness;
		
	};


}

#endif // !OFEC_CSIWDN_Solution_H


