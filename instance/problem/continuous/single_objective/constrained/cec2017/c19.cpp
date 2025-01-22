#include "C19.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {

		void C19::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-50., 50.);
			setInitialDomain(-50., 50.);
			m_num_cons = 2;
			m_constraint.resize(2);
			for (auto& i : m_constraint)
				i = Constraint::kInequality;

			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C19::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];
			
			size_t i;
			obj[0] = 0.;

			for (i = 0; i<m_num_vars; ++i)
			{
				obj[0] += (pow(fabs(x[i]), 0.5) + 2 * sin(pow(x[i], 3)));
			}
			obj[0] += m_bias;

			// evaluate constraint value

			for (auto &i : con)
				i = 0.;

			for (i = 0; i<m_num_vars - 1; ++i)
			{
				con[0] += -10 * exp(-0.2*sqrt(x[i] * x[i] + x[i + 1] * x[i + 1]));
			}
			con[0] += (m_num_vars - 1) * 10 / exp(-5.0);
			if (con[0] <= 0) con[0] = 0;

			for (i = 0; i<m_num_vars; ++i)
			{
				con[1] += sin(2 * x[i])*sin(2 * x[i]);
			}
			con[1] -= 0.5 * m_num_vars;
			if (con[1] <= 0) con[1] = 0;
			
		}
	}
}