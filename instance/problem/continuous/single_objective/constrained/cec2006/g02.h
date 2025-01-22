/********* Begin Register Information **********
{"name": "COP_CEC2006_G02","keyword": "CEC2006_COP_G02","problem tags": [ "ConOP", "COP", "SOP" ]}
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


#ifndef OFEC_G02_H
#define OFEC_G02_H

#include "../cop_base.h"

namespace ofec {
	namespace CEC2006 {
		class G02 final : public CopBase
		{
		protected:
			void initialize_() override;
			void evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) override;
		};
	}
	using CEC2006_COP_G02 = CEC2006::G02;
}
#endif // ! OFEC_G02_H

