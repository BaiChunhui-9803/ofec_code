#include "C03.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {
		
		void C03::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
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

		void C03::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			size_t i, j;
			obj[0] = 0.;
			for (i = 1; i < m_num_vars + 1; ++i)
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

			con[0] = 0;
			for (i = 0; i < m_num_vars; ++i)
			{
				con[0] += (x[i] * x[i] - 5000 * cos(0.1*OFEC_PI*x[i]) - 4000);
			}
			if (con[0] <= 0) con[0] = 0;

			con[1] = 0;
			for (i = 0; i < m_num_vars; ++i)
			{
				con[1] -= x[i] * sin(0.1*OFEC_PI*x[i]);
			}
			if (fabs(con[1]) - 1e-4 <= 0) con[1] = 0;
		}
	}
}