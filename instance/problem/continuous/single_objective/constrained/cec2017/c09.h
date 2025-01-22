/********* Begin Register Information **********
{
	"name": "COP_CEC2017_C09",
	"keyword": "CEC2017_COP_C09",
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


#ifndef OFEC_C09_H
#define OFEC_C09_H

#include "../cop_base.h"

namespace ofec {
	namespace CEC2017 {
		class C09 final: public CopBase
		{
		public:
			void initialize_();
		protected:
			void evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) override;
		private:

		};
	}
	using CEC2017_COP_C09 = CEC2017::C09;
}
#endif // ! OFEC_C09_H





