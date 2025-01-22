/********* Begin Register Information **********
{
	"name": "NBN-COM-ALG",
	"identifier": "NBN_COM_ALG",
	"problem tags": [ "ComOP", "GOP", "SOP", "MMOP", "TSP" ]
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

#ifndef OFEC_NBN_COM_ALG_H
#define OFEC_NBN_COM_ALG_H

#include "../../../../core/algorithm/algorithm.h"
#include "nbn_alg_com_gl_pop.h"
namespace ofec {
	class NBN_COM_ALG : public Algorithm {
	protected:
		std::unique_ptr<PopGL_NBN_COM_ALG>  m_pop;
		size_t m_pop_size;

		void run_() override;
		void initialize_() override;
#ifdef OFEC_DEMO
		void updateBuffer();
#endif

	public:
		//void record() override;
	};
}

#endif // !OFEC_CMAES_H

