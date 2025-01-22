/********* Begin Register Information **********
{
	"name": "EOP_CEC2015_F15",
	"identifier": "CEC2015_EOP_F15",
	"problem tags": [ "EOP", "SOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

#ifndef OFEC_F15_EXPENSIVE_COMPOSITION3_H
#define OFEC_F15_EXPENSIVE_COMPOSITION3_H

#include "../../global/cec2015/composition.h"

namespace ofec {
	namespace cec2015 {
		class F15_expensive_composition3 : public composition_2015
		{
		public:
			F15_expensive_composition3(const ParameterMap &v);
			F15_expensive_composition3(const std::string &name, size_t size_var, size_t size_obj);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			void setFunction();

		private:

		};
	}
	using CEC2015_EOP_F15 = CEC2015::F15_expensive_composition3;
}
#endif // !OFEC_F15_EXPENSIVE_COMPOSITION3_H


