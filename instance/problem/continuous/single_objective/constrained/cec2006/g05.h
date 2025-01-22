/********* Begin Register Information **********
{
	"name": "COP_CEC2006_G05",
	"keyword": "CEC2006_COP_G05",
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


#ifndef OFEC_G05_H
#define OFEC_G05_H

#include "../cop_base.h"

namespace ofec {
	namespace CEC2006 {
		class G05 final : public CopBase
		{
		public:
			void initialize_() override;
		protected:
			void evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) override;
		private:

		};
	}
	using CEC2006_COP_G05 = CEC2006::G05;
}
#endif // ! OFEC_G05_H

