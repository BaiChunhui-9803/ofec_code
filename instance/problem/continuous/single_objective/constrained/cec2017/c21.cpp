#include "C21.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {

		void C21::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			m_num_cons = 2;
			m_constraint.resize(2);
			m_constraint[0] = Constraint::kInequality;
			m_constraint[1] = Constraint::kEquality;

			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path
			loadRotation("/instance/problem/continuous/constrained/CEC2017/data/");
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C21::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			size_t i;

			obj[0] = 0.;
			for (i = 0; i < m_num_vars; ++i)
			{
				obj[0] += (x[i] * x[i] - 10.0 * cos(2.0 * OFEC_PI * x[i]) + 10.0);
			}

			obj[0] += m_bias;

			// evaluate constraint value

			std::vector<Real> x_ro(x, x + m_num_vars);
			rotate(x_ro.data());

			for (auto& i : con)
				i = 0.;

			for (i = 0; i < m_num_vars; ++i)
			{
				con[0] += fabs(x_ro[i]);
			}
			con[0] = 4 - con[0];
			if (con[0] <= 0) con[0] = 0;

			for (i = 0; i < m_num_vars; ++i)
			{
				con[1] += x_ro[i] * x_ro[i];
			}
			con[1] -= 4;
			if (fabs(con[1]) - 1e-4 <= 0) con[1] = 0;

		}
	}
}