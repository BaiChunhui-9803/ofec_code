/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F02",
	"identifier": "GOP_CEC2005_F02",
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


#ifndef OFEC_F2_SHIFTED_SCHWEFEL_1_2_H
#define OFEC_F2_SHIFTED_SCHWEFEL_1_2_H

#include "../classical/schwefel_1_2.h"

namespace ofec {
	namespace cec2005 {
		class ShiftedSchwefel_1_2 : public Schwefel_1_2 {
			OFEC_CONCRETE_INSTANCE(ShiftedSchwefel_1_2)
		protected:
			void addInputParameters() {}
			void initialize_(Environment *env) override;
			void updateOptima(Environment *env) override;
		};
	}
	using GOP_CEC2005_F02 = cec2005::ShiftedSchwefel_1_2;
}
#endif // ! OFEC_F2_SHIFTED_SCHWEFEL_1_2_H
