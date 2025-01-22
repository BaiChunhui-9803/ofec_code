/********* Begin Register Information **********
{
	"name": "GOP_CEC2015_F03",
	"identifier": "CEC2015_GOP_F03",
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
#ifndef OFEC_F3_SR_ACKLEY_H
#define OFEC_F3_SR_ACKLEY_H

#include "../classical/ackley.h"

namespace ofec {
	namespace cec2015 {
		class F3_SR_ackley final : public Ackley {
		protected:
			void initialize_(Environment *env) override;
		};
	}
	using CEC2015_GOP_F03 = cec2015::F3_SR_ackley;
}

#endif // !OFEC_F3_SR_ACKLEY_H



