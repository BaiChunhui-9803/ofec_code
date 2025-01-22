/********* Begin Register Information **********
{
	"name": "Classic_center_peak",
	"identifier": "CenterPeak",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper: Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
* 
//Ursem F4 in R. K. Ursem, ��Multinational evolutionary algorithms,�� in Proc. Congr.
//Evol. Comput. (CEC), vol. 3. 1999, pp. 1633�C1640.
//Ursem F4
*******************************************************************************************/

#ifndef OFEC_CENTER_PEAK_H
#define OFEC_CENTER_PEAK_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	class CenterPeak : public Continuous, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(CenterPeak)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
	};
}
#endif // !OFEC_CENTER_PEAK_H