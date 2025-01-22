#include "C14.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {

		void C14::initialize_() {
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
		void C14::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];
			

			size_t i;
			Real num1=0., num2=0.;
			obj[0] = 0.;
			for (i = 0; i<m_num_vars; ++i)
			{
				num1 += x[i] * x[i];
			}
			num1 = sqrt(num1 / m_num_vars);
			for (i = 0; i<m_num_vars; ++i)
			{
				num2 += cos(2 * OFEC_PI*x[i]);
			}
			num2 /= m_num_vars;
			obj[0] = -20 * exp(-0.2*num1) + 20 - exp(num2) + exp(1.0);
			obj[0] += m_bias;

			
			// evaluate constraint value

			for (auto &i : con)
				i = 0.;

			for (i = 1; i < m_num_vars; ++i)
			{
				con[0] += pow(x[i], 2);
			}
			con[0] += 1 - fabs(x[0]);
			if (con[0] <= 0) con[0] = 0;

			for (i = 0; i < m_num_vars; ++i)
			{
				con[1] += x[i] * x[i];
			}
			con[1] -= 4;
			if (fabs(con[1]) - 1e-4 <= 0) con[1] = 0;
			
		}
	}
}