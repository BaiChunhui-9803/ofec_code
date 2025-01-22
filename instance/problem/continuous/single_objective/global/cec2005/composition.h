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

/*------------------------------------------------------------------------------------------------------------
Suganthan, P. N., Hansen, N., Liang, J. J., Deb, K., Chen, Y. P., Auger, A., & Tiwari, S. (2005).
Problem definitions and evaluation criteria for the CEC 2005 special session on Real-parameter optimization.
KanGAL report, 2005005, 2005.
------------------------------------------------------------------------------------------------------------*/

#ifndef OFEC_CEC2005_COMPOSITION_H
#define OFEC_CEC2005_COMPOSITION_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../function.h"

namespace ofec::cec2005 {
	class Composition : virtual public Continuous {
		OFEC_ABSTRACT_INSTANCE(Composition)
	public:
		size_t numFunctions();
		Function* getFunction(size_t num);

	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;

		void computeFmax();
		virtual void setFunction(Environment *env) = 0;
		virtual void setWeight(std::vector<Real> &w, const std::vector<Real> &x) const;
		virtual bool loadTranslation(const std::string &path);
		virtual void setTranslation();
		virtual bool loadRotation(const std::string &path);
		virtual void setRotation();

	protected:
		size_t m_num_function;                // number of basic functions, for hybrid functions
		std::vector<std::shared_ptr<Function>> m_function;    // the functions
		std::vector<std::shared_ptr<ParameterMap>> m_param_fun;		  // param of functions
		std::vector<Real> m_height;
		std::vector<Real> m_fmax;
		Real m_height_normalize_severity;   // constant number for noralizing all basic function with similar height
		std::vector<Real> m_converge_severity;     // severity of converge range for each function
		std::vector<Real> m_stretch_severity;      // severity of stretching original function, greater than 1 for stretch
		std::string m_index_number;			// index number in the benchmark suite
	};
}
#endif // ! OFEC_CEC2005_COMPOSITION_H
