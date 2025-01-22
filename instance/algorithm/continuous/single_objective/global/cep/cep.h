/********* Begin Register Information **********
{
	"name": "CEP",
	"identifier": "CEP",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*-------------------------------------------------------------------------------
*
* Created on 15th Aug, 2019 by Junchen Wang
*
*********************************************************************************/

#ifndef OFEC_CEP_H
#define OFEC_CEP_H

#include "../../../../template/classic/evolutionary_programming/population.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class CEP : virtual public Algorithm {   
		OFEC_CONCRETE_INSTANCE(CEP)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

	protected:
		Real m_tau, m_tau_prime;
		size_t m_q, m_pop_size;
		PopEP<> m_pop;
	};
}

#endif // !OFEC_CEP_H
