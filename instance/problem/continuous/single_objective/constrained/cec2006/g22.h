/********* Begin Register Information **********
{
	"name": "COP_CEC2006_G22",
	"keyword": "CEC2006_COP_G22",
	"problem tags": [ "ConOP", "COP", "SOP" ]
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


#ifndef OFEC_G22_H
#define OFEC_G22_H

#include "../cop_base.h"

namespace ofec {
	namespace CEC2006 {
		class G22 final : public CopBase
		{
		public:
			void initialize_() override;
		protected:
			void evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) override;
		private:

		};
	}
	using CEC2006_COP_G22 = CEC2006::G22;
}
#endif // ! OFEC_G22_H


