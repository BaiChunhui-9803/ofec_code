/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F14",
	"identifier": "GOP_CEC2005_F14",
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


#ifndef OFEC_F14_SHIFTED_ROTATED_EXPANDED_SCAFFER_F6_H
#define OFEC_F14_SHIFTED_ROTATED_EXPANDED_SCAFFER_F6_H

#include "../classical/scaffer_F6.h"

namespace ofec {
	namespace cec2005 {
		class ShiftedRotatedExpandedScafferF6 : public ScafferF6 {
			OFEC_CONCRETE_INSTANCE(ShiftedRotatedExpandedScafferF6)
		protected:
			void addInputParameters() {}
			void initialize_(Environment *env) override;
			void updateOptima(Environment *env) override;
		};
	}
	using GOP_CEC2005_F14 = cec2005::ShiftedRotatedExpandedScafferF6;
}
#endif // ! OFEC_F14_SHIFTED_ROTATED_EXPANDED_SCAFFER_F6_H


