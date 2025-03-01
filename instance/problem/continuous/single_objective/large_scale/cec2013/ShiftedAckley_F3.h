/********* Begin Register Information **********
{
	"name": "LSOP_CEC2013_F03",
	"identifier": "CEC2013_LSOP_F03",
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
// Created: 8 August 2017
// Last modified:
#ifndef _SHIFTEDACKLEY_F3_H
#define _SHIFTEDACKLEY_F3_H

#include "function_CEC2013.h"
namespace ofec {
	namespace CEC2013 {
		class ShiftedAckley_F3 final : public function_CEC2013 {
		public:
			ShiftedAckley_F3(const ParameterMap &v);
			ShiftedAckley_F3(const std::string &name, size_t size_var, size_t size_obj);
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			~ShiftedAckley_F3();

			void initialize();
		};
	}
	using CEC2013_LSOP_F03 = CEC2013::ShiftedAckley_F3;
}
#endif
