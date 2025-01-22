#include "C02.h"
#include "../../../../../core/instance_manager.h"


namespace ofec {
	namespace CEC2017 {

		void C02::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			m_num_cons = 1;
			m_constraint.resize(1);
			m_constraint[0] = Constraint::kInequality;

			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/COP_CEC2017_C02_Shift_2Dim");  //data path
			loadRotation("/instance/problem/continuous/constrained/CEC2017/data/COP_CEC2017_C02_RotM_2Dim");

			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}

		void C02::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			size_t i, j;
			obj[0] = 0.;
			for (i = 1; i<m_num_vars+1 ; ++i)
			{
				Real t = 0.0;
				for (j = 0; j<i; j++)
				{
					t += x[j];
				}
				obj[0] += t*t;
			}
			obj[0] += m_bias;

			// evaluate constraint value

			std::vector<Real> x_(x, x + m_num_vars);
			rotate(x_.data());
			con[0] = 0;
			for (i = 0; i < m_num_vars; ++i)
			{
				con[0] += (x_[i] * x_[i] - 5000 * cos(0.1*OFEC_PI*x_[i]) - 4000);
			}
			if (con[0] <= 0) con[0] = 0;
		}
	}
}