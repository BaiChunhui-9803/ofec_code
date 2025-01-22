/********* Begin Register Information **********
{
	"name": "LSOP_CEC2013_F14",
	"identifier": "CEC2013_LSOP_F14",
	"problem tags": [ "LSOP", "SOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Li Zhou
* Email: 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

#ifndef CONFLICTINGOS_SHIFTEDSCHWEFEL_F14_H
#define CONFLICTINGOS_SHIFTEDSCHWEFEL_F14_H

#include "function_CEC2013.h"
namespace ofec {
	namespace CEC2013 {
		class ConflictingOS_ShiftedSchwefel_F14 :public function_CEC2013 {
		public:
			ConflictingOS_ShiftedSchwefel_F14(const ParameterMap &v);
			ConflictingOS_ShiftedSchwefel_F14(const std::string &name, size_t size_var, size_t size_obj);
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			~ConflictingOS_ShiftedSchwefel_F14();

			void initialize();
		};
	}
	using CEC2013_LSOP_F14 = CEC2013::ConflictingOS_ShiftedSchwefel_F14;
}
#endif





