/********* Begin Register Information **********
{
	"name": "GOP_CEC2015_F13",
	"identifier": "CEC2015_GOP_F13",
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
*************************************************************************
*  Paper : Problem Definitions and Evaluation Criteria for the CEC2015
*  Competition on Learning-based Real-Parameter Single Objective
*  Optimization.
************************************************************************/

#ifndef OFEC_F13_GLOBAL_COMPOSITION5_H
#define OFEC_F13_GLOBAL_COMPOSITION5_H

#include "composition.h"
#include "hybrid.h"

namespace ofec {
	namespace cec2015 {
		class F13_global_composition5 final : public composition_2015
		{
		public:
			F13_global_composition5(const ParameterMap &v);
			F13_global_composition5(const std::string &name, size_t size_var, size_t size_obj);
			~F13_global_composition5();
			Hybrid* get_hybrid(size_t num);
			void initialize();
		protected:
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			void setFunction();
			bool load_hybrid_translation(const std::string &path);
			void set_hybrid_translation();
			void set_weight(std::vector<Real>& weight, const std::vector<Real>&x);
		protected:
			size_t m_num_hybrid;
			std::vector<Hybrid*> m_hybrid;
		};
	}
	using CEC2015_GOP_F13 = cec2015::F13_global_composition5;
}
#endif // !OFEC_F13_GLOBAL_COMPOSITION5_H




