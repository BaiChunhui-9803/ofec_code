/********* Begin Register Information **********
{
	"name": "Classic_FIBA",
	"identifier": "FIBA",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (ofec)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of ofec. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper: A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/

#ifndef OFEC_FIBA_H
#define OFEC_FIBA_H

#include "../../function.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {	
	class FIBA : public Function, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(FIBA)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateOriginalObj(Real *x, std::vector<Real>& obj) const override;
		void setPara();

		Real m_kappa, m_chi;
		int m_case;

	public:
		void setCase(int c);
	};
	
}
#endif // !OFEC_FIBA_H