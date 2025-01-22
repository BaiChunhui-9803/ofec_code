/********* Begin Register Information **********
{
	"name": "GOP_CEC2015_F10",
	"identifier": "CEC2015_GOP_F10",
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

#ifndef OFEC_F10_GLOBAL_COMPOSITION2_H
#define OFEC_F10_GLOBAL_COMPOSITION2_H

#include "composition.h"
#include "hybrid.h"

namespace ofec {
	namespace cec2015 {
		class F10_composition2 final : public composition_2015 {
		public:
			~F10_composition2();
			Hybrid* getHybrid(size_t num);	
		protected:
			void initialize_(Environment *env) override;
			void evaluateObjective(Real *x, std::vector<Real>& obj) override;
			void setFunction() override;
			bool loadTranslation(const std::string &path) override;
			void setTranslation() override;
			void setWeight(std::vector<Real>& weight, const std::vector<Real>&x) override;
		protected:
			std::vector<Hybrid*> m_hybrid;
			size_t m_num_hybrid;
			std::vector<Real> m_hy_bias;
		};
	}
	using CEC2015_GOP_F10 = cec2015::F10_composition2;
}
#endif // !OFEC_F10_GLOBAL_COMPOSITION2_H


