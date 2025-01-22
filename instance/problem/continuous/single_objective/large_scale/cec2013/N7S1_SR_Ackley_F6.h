/********* Begin Register Information **********
{
	"name": "LSOP_CEC2013_F06",
	"identifier": "CEC2013_LSOP_F06",
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

#ifndef N7S1_SR_ACKLEY_F6_H
#define N7S1_SR_ACKLEY_F6_H

#include "function_CEC2013.h"
namespace ofec {
	namespace CEC2013 {
		class N7S1_SR_Ackley_F6 final:public function_CEC2013 {
		public:
			N7S1_SR_Ackley_F6(const ParameterMap &v);
			N7S1_SR_Ackley_F6(const std::string &name, size_t size_var, size_t size_obj);
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			~N7S1_SR_Ackley_F6();

			void initialize();
		};
	}
	using CEC2013_LSOP_F06 = CEC2013::N7S1_SR_Ackley_F6;
}
#endif

