/********* Begin Register Information **********
{
	"name": "COP_CEC2017_C05",
	"keyword": "CEC2017_COP_C05",
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


#ifndef OFEC_C05_H
#define OFEC_C05_H

#include "cop_base.h"

namespace ofec {
	namespace CEC2017 {
		class C05 final: public CopBase
		{
		public:
			void initialize_();
		protected:
			void evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) override;
			bool load_rotation_C05(const std::string &path);
			void load_rotation_C05_(const std::string &path);
			void set_rotation_C05();
			void rotate(Real *x, size_t num);
		private:
			Matrix m_mat1, m_mat2;
		};
	}
	using CEC2017_COP_C05 = CEC2017::C05;
}
#endif // ! OFEC_C05_H


