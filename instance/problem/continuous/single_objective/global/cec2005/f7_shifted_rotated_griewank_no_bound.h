/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F07",
	"identifier": "GOP_CEC2005_F07",
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


#ifndef OFEC_F7_SHIFTED_ROTATED_GRIEWANK_NO_BOUND_H
#define OFEC_F7_SHIFTED_ROTATED_GRIEWANK_NO_BOUND_H

#include "../classical/griewank.h"
#include "../../../init_pop_bounded.h"

namespace ofec {
	namespace cec2005 {
		class ShiftedRotatedGriewankNoBound : public InitPopBounded, public Griewank {
			OFEC_CONCRETE_INSTANCE(ShiftedRotatedGriewankNoBound)
		protected:
			void addInputParameters() {}
			void initialize_(Environment *env) override;
			void updateOptima(Environment *env) override;
		};
	}
	using GOP_CEC2005_F07 = cec2005::ShiftedRotatedGriewankNoBound;
}
#endif // ! OFEC_F7_SHIFTED_ROTATED_GRIEWANK_NO_BOUND_H


