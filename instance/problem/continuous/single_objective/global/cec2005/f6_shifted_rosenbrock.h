/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F06",
	"identifier": "GOP_CEC2005_F06",
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
*************************************************************************/


#ifndef OFEC_F6_SHIFTED_ROSENBROCK_H
#define OFEC_F6_SHIFTED_ROSENBROCK_H

#include "../classical/rosenbrock.h"

namespace ofec {
	namespace cec2005 {
		class ShiftedRosenbrock : public Rosenbrock {
			OFEC_CONCRETE_INSTANCE(ShiftedRosenbrock)
		protected:
			void addInputParameters() {}
			void initialize_(Environment *env) override;
			void updateOptima(Environment *env) override;
		};
	}
	using GOP_CEC2005_F06 = cec2005::ShiftedRosenbrock;
}
#endif // ! OFEC_F6_SHIFTED_ROSENBROCK_H

