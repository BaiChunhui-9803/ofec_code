/********* Begin Register Information **********
{
	"name": "SaDE",
	"identifier": "SaDE",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

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
// updated Mar 28, 2018 by Li Zhou

/*
Paper: A. K. Qin, V. L. Huang and P. N. Suganthan, ��Differential evolution with 
strategy adaptation for global numerical optimization,�� IEEE Transactions on 
Evolutionary Computation, 13(2): 398- 417, 2009. 
*/
#ifndef OFEC_SADE_H
#define OFEC_SADE_H

#include "sade_pop.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class SaDE : virtual public Algorithm	{
		OFEC_CONCRETE_INSTANCE(SaDE)
	protected:
		std::unique_ptr<PopSaDE> m_pop;
		size_t m_pop_size;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

	public:		
		const std::vector<Real>& ratioStrategy() const { return m_pop->ratioStrategy(); }
	};
}
#endif // OFEC_SADE_H
