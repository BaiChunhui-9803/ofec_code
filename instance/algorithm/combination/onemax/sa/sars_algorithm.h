/********* Begin Register Information **********
{
	"name": "SARS-OneMax",
	"identifier": "SARS_OneMax",
	"problem tags": [ "ComOP", "GOP", "SOP", "MMOP", "OneMax" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*-------------------------------------------------------------------------------
*
*
*-------------------------------------------------------------------------------
*
*
*************************************************************************/


#ifndef OFEC_SARS_ONEMAX_H
#define OFEC_SARS_ONEMAX_H
#include "sars.hpp"
#include "../../../../../core/algorithm/algorithm.h"


namespace ofec {
	class SARS_OneMax : public Algorithm {

	protected:
		int m_nvar = 0;
		double p = 0;
		double Tend = 0;
		double epsilon = 0;
		double innerBudgetParam = 0;

		unsigned long long int step = 1;
		unsigned long long int innerBudget = 512;
		unsigned long long int stepMain = 1;


		double Tstart = 0;
		//T_from_DeltaE_and_P(max(1.0, n / 4.0), 0.1);
		std::unique_ptr<SolBase> m_cur_sol;
		std::unique_ptr<SolBase> m_best_sol;
		std::unique_ptr<SolBase> m_new_sol;

		void evolve();
		
	public:

		void run_() override;
		void initialize_() override;


	};

}

#endif // !OFEC_SARS_ONEMAX_H
