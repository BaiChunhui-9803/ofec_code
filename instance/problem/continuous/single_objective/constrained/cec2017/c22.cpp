#include "C22.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {
	
		void C22::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			m_num_cons = 3;
			m_constraint.resize(3);
			for (auto& i : m_constraint)
				i = Constraint::kInequality;

			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path

			loadRotation("/instance/problem/continuous/constrained/CEC2017/data/");
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C22::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			std::vector<Real> x_ori(x, x + m_num_vars);
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];
			rotate(x);

			size_t i;

			obj[0] = 0.;
			for (i = 0; i<m_num_vars - 1; ++i)
			{
				obj[0] += (100 * (x[i] * x[i] - x_ori[i + 1])*(x[i] * x[i] - x_ori[i + 1]) + (x[i] - 1)*(x[i] - 1));
			}
			obj[0] += m_bias;

			// evaluate constraint value

			for (auto &i : con)
				i = 0.;

			for (i = 0; i < m_num_vars; ++i)
			{
				con[0] += x[i] * x[i] - 10 * cos(2 * OFEC_PI*x[i]) + 10;
			}
			con[0] -= 100;
			if (con[0] <= 0) con[0] = 0;

			for (i = 0; i < m_num_vars; ++i)
			{
				con[1] += x[i];
			}
			con[1] -= 2 * m_num_vars;
			if (con[1] <= 0) con[1] = 0;
			
			for (i = 0; i < m_num_vars; ++i)
			{
				con[2] += x[i];
			}
			con[2] = 5 - con[2];
			if (con[2] <= 0) con[2] = 0;
		}
	}
}