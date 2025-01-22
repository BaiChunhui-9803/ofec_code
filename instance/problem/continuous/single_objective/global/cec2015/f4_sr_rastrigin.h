/********* Begin Register Information **********
{
	"name": "GOP_CEC2015_F04",
	"identifier": "CEC2015_GOP_F04",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

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
*************************************************************************
*  Paper : Problem Definitions and Evaluation Criteria for the CEC2015
*  Competition on Learning-based Real-Parameter Single Objective
*  Optimization.
************************************************************************/
#ifndef OFEC_F4_SR_RASTRIGIN_H
#define OFEC_F4_SR_RASTRIGIN_H

#include "../classical/rastrigin.h"

namespace ofec {
	namespace cec2015 {
		class F4_SR_rastrigin final : public Rastrigin {
			void initialize_(Environment *env) override;
		};
	}
	using CEC2015_GOP_F04 = cec2015::F4_SR_rastrigin;
}

#endif // !OFEC_F4_SR_RASTRIGIN_H




