/********* Begin Register Information **********
{
	"name": "LSOP_CEC2013_F01",
	"identifier": "CEC2013_LSOP_F01",
	"problem tags": [ "LSOP", "SOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

#ifndef OFEC_SHIFTEDELLIPTIC_F1_H
#define OFEC_SHIFTEDELLIPTIC_F1_H

#include "function_CEC2013.h"

namespace ofec {
	namespace CEC2013 {
		class ShiftedElliptic_F1 final :public function_CEC2013 {
		public:
			ShiftedElliptic_F1(const ParameterMap &v);
			ShiftedElliptic_F1(const std::string &name, size_t size_var, size_t size_obj);
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			~ShiftedElliptic_F1();

			void initialize();
		};
	}
	using CEC2013_LSOP_F01 = CEC2013::ShiftedElliptic_F1;
}
#endif // ! OFEC_SHIFTEDELLIPTIC_F1_H