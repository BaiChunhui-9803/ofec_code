#include "C13.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {

		void C13::initialize_() {
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
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C13::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			size_t i;
			obj[0] = 0.;
			for (i = 0; i < m_num_vars - 1; ++i)
			{
				obj[0] += (100 * (x[i] * x[i] - x[i + 1])*(x[i] * x[i] - x[i + 1]) + (x[i] - 1)*(x[i] - 1));
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
			for (i = 0; i < m_num_vars; ++i)
			{
				con[1] += x[i];
			}
			con[1] -= 2 * m_num_vars;
			for (i = 0; i < m_num_vars; ++i)
			{
				con[2] += x[i];
			}
			con[2] = 5 - con[2];

			for (auto &i : con)
				if (i <= 0) i = 0;
		}
	}
}