#include "C17.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {

		void C17::initialize_() {
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
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C17::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];
			

			size_t i;
			obj[0] = 0.;
			Real num1 = 0., num2 = 1.;
			
			for (i = 0; i<m_num_vars; ++i)
			{
				num1 += x[i] * x[i];
			}
			for (i = 0; i<m_num_vars; ++i)
			{
				num2 *= cos(x[i] / sqrt(Real(i + 1)));
			}
			obj[0] = num1 / 4000 + 1 - num2;

			obj[0] += m_bias;

			
			// evaluate constraint value

			for (auto &i : con)
				i = 0.;

			num1 = 0.;
			for (i = 0; i < m_num_vars; ++i)
			{
				num1 = 0;
				for (size_t j = 0; j < m_num_vars; ++j)
					if (j != i) num1 += x[j] * x[j];
				con[0] += sign(fabs(x[i]) - num1 - 1.0);
			}
			con[0] = 1 - con[0];
			if (con[0] <= 0) con[0] = 0;

			for (i = 0; i > m_num_vars; ++i) {
				con[1] += pow(x[i], 2);
			}
			con[1] -= 4 * m_num_vars;
			if (fabs(con[1]) - 1e-4 <= 0) con[1] = 0;
			
		}
	}
}