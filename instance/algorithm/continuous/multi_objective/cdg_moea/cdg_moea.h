/********* Begin Register Information **********
{
	"name": "CDG-MOEA",
	"identifier": "CDGMOEA",
	"problem tags": [ "ConOP", "MOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Xiaofang Wu
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************
*  See the details of CDG-MOEA in the following paper
*  X.Cai, Z.Mei, Z.Fan, et al. A Constrained Decomposition Approach with Grids for Evolutionary 
*  Multiobjective Optimization. IEEE Transactions on Evolutionary Computation, 2017:1-1.
*************************************************************************/
// Created: 4 May 2019
// Last modified: 24 Aug 2019 by Xiaofang Wu (email:wuxiaofang@cug.edu.cn)

#ifndef OFEC_CDG_MOEA_H
#define OFEC_CDG_MOEA_H

#include "../../../template/classic/de/individual.h"
#include "../../../template/classic/de/population.h"
#include "../../../../../core/algorithm/algorithm.h"
#include "../moea_de_pop.h"
#include <algorithm>

namespace ofec {
	/****************************** Solution *************************/
	class IndCDGMOEA : public IndDE {
		std::vector<int> m_Gk;
		std::vector<size_t> m_GN;
		std::vector<int> m_R;
		double m_meta = 20;
	public:
		IndCDGMOEA(size_t number_objectives, size_t num_cons, size_t num_vars) :
			IndDE(number_objectives, num_cons, num_vars),
			m_Gk(number_objectives),
			m_R(number_objectives, 0) {}
						
		std::vector<int>& Gk() { return m_Gk; }
		const std::vector<int>& Gk()const { return m_Gk; }
		std::vector<size_t>& GN() { return m_GN; }
		const std::vector<size_t>& GN() const { return m_GN; }
		std::vector<int>& R() { return m_R; }
		const std::vector<int>& R() const { return m_R; }
	};

	/***************************** population *************************/
	class PopCDGMOEA : public PopMODE<IndCDGMOEA> {
	protected:
		std::vector<double> m_ideal;
		std::vector<double> m_nadir;
		std::vector<std::pair<Real, Real>> m_grid_min_max;	//first is m_grid_min,and second is m_grid_max
		std::vector<Real> m_grid_distance;
		int m_grid_div;	//division numbers K=180 for bi-objective, K=30 for tri-objective
		int m_T;	//neighborhood distance set 5 for bi-objective and 1 for tri-objective
		double m_delta = 0.9;
		std::vector<std::vector<std::vector<int>>> m_S;
		std::vector<std::vector<int>> m_R;
		std::vector<IndCDGMOEA> m_offspring;	//store P and Q			

	public:
		PopCDGMOEA(size_t size_pop,Problem *pro);
		void initialize(Problem *pro, Random *rnd) override;
		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;

		std::vector<double>& ideal() { return m_ideal; }
		const std::vector<double>& ideal() const { return m_ideal; }
		std::vector<double>& nadir() { return m_nadir; }
		const std::vector<double>& nadir() const { return m_nadir; }
		std::vector<std::vector<std::vector<int>>>& S() { return m_S; }
		const std::vector<std::vector<std::vector<int>>>& S() const  { return m_S; }
		std::vector<std::vector<int>>& R() { return m_R; }
		const std::vector<std::vector<int>>& R() const { return m_R; }
		void updateIdealPoint(Problem *pro);
		void updateIdealPoint(const std::vector<IndCDGMOEA> &offspring,Problem *pro);
		void updateNadirPoint(Problem *pro);
		void updateNadirPoint(const std::vector<IndCDGMOEA> &offspring,Problem *pro);
		
		void gridConstruct(Problem *pro);	// construct grid environment
		void assign_GN(Problem *pro);
		void gridConstruct_assignGN_P_reserve(std::vector<IndCDGMOEA> &offspring, std::vector<int>& P_reserve, int size,Problem *pro);
		void assign_S(std::vector<IndCDGMOEA> &offspring, std::vector<int> P_reserve,Problem *pro);
		void RBS(std::vector<IndCDGMOEA> &offspring, std::vector<int> P_reserve,Problem *pro);

		void selectRandom(int number, std::vector<int>& result, Random *rnd);
		void nondominatedSorting(std::vector<IndCDGMOEA> &offspring,Problem *pro);
			
		//void record();
		//void record_x();
		//void record_x_offspring();
		//void record_f();
		//void record_f_offspring();
		int evolveMO(Problem *pro,Algorithm *alg, Random *rnd);
		void evalEens(Problem *pro, Random *rnd);

	private:
		void sortl(std::vector<std::pair<double, int>>& v);
		void sortR(const std::vector<int>& v1, std::vector<int>& v2);
		void sortL(std::vector<std::vector<int>>& R, std::vector<int>& number,Problem *pro);			
	};		

	/******************************CDGMOEA****************************/
	class CDGMOEA : public Algorithm {
	public:
		void initialize_() override;
		void record() override;
		void initPop();
#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;
		std::unique_ptr<PopCDGMOEA> m_pop;
		size_t m_pop_size;
	};

}


#endif
