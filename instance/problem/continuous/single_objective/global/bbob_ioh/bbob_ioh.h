/********* Begin Register Information **********
[
	{ "name":"BBOB_IOH", "identifier":"BBOB_IOH", "problem tags":["ConOP","GOP","SOP"] }
]
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li and Yiya Diao
* Email: changhe.lw@gmail.com, 744862398@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// source code refered to bbob.v15.02,https://github.com/IOHprofiler/IOHexperimenter
#ifndef OFEC_BBOB_IOH_H
#define OFEC_BBOB_IOH_H

#include "../../../../../../core/problem/solution.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_gop.h"
#include "../../../../../../utility/IOHexperimenter/include/ioh/problem/bbob.hpp"
#include <string>
#include <algorithm>

namespace ofec {

#define GET_BBOB_IOH(pro) dynamic_cast<BBOB_IOH*>(pro)


	class BBOB_IOH final : public Continuous, public MetricsGOP {
	protected:

		//ioh::common::Factory< ioh::problem::BBOB, int, int>  m_problem_factory;
		std::shared_ptr<ioh::problem::BBOB>	m_ptr_problem;

		int m_func_id = 1;
		int m_ins_id = 1;
	protected:
		void initialize_() override {

			Continuous::initialize_();
			auto& v = *m_param;
			resizeVariable(v.get<int>("number of variables"));
			resizeObjective(1);
			m_optimize_mode[0] = OptimizeMode::kMinimize;

			const auto& problem_factory = ioh::problem::ProblemRegistry<ioh::problem::BBOB>::instance();
			m_ptr_problem = problem_factory.create(m_func_id, m_ins_id, m_number_variables);
			
			const auto& bounds = m_ptr_problem->bounds();
			m_domain.resize(m_number_variables);
			for (size_t i = 0; i < m_domain.size(); ++i)
				m_domain.setRange(bounds.lb[i], bounds.ub[i], i);
			m_domain_update = true;
			m_domain_volume_update = true;
		
			m_optima.reset(new Optima<>());
			dynamic_cast<Optima<VariableVector<>>&>(*m_optima).appendVar(VariableVector<>(m_ptr_problem->optimum().x));
			m_optima->setVariableGiven(true);
			m_optima->appendObj(m_ptr_problem->optimum().y);
			m_optima->setObjectiveGiven(true);

		}
	

	private:

	public:


		void setFuncId(int func_id) {
			m_func_id = func_id;
			
		}
		void setInsId(int ins_id) {
			m_ins_id = ins_id;
		}
		virtual void evaluate_(SolutionBase& s, bool effective) override{
			VariableVector<Real>& x = dynamic_cast<Solution<>&>(s).variable();
			auto& obj = dynamic_cast<Solution<> &>(s).objective();
			auto& con = dynamic_cast<Solution<> &>(s).constraint();

			std::vector<double> x_(x.begin(), x.end()); //for parallel running
			obj.resize(m_number_objectives);
			obj.front()= (*m_ptr_problem)(x_);
		}
		void test() {
			const auto& problem_factory = ioh::problem::ProblemRegistry<ioh::problem::BBOB>::instance();
			int func_id = 1;
			int ins_id = 1;
			std::vector<double> x = { -4.999921736,-3.684622119,2.556053222,-0.4134986808,0.3276723741 };
			double f = 138.7457513;

			auto problem = problem_factory.create(func_id, ins_id, 5);
			auto y = (*problem)(x);
	}
	};
}
#endif // !OFEC_BBOB_H
