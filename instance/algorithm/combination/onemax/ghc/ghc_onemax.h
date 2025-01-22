/********* Begin Register Information **********
{
	"name": "GHC-OneMax",
	"identifier": "GHC_OneMax",
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


#ifndef OFEC_GHC_ONEMAX_H
#define OFEC_GHC_ONEMAX_H



#include "../../../../../core/algorithm/algorithm.h"


namespace ofec {
	class GHC_OneMax : public Algorithm {

	protected:
		std::unique_ptr<SolBase> m_cur_sol;
		std::unique_ptr<SolBase> m_best_sol;
		int m_flit_idx = 0;
		int m_dim = 0;
		
		void evolve();
	
	public:

		void run_() override;
		void initialize_() override;
	};

}

#endif // !OFEC_SARS_ONEMAX_H
