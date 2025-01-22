#include "C25.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {
	
		void C25::initialize_() {
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
		void C25::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];
			rotate(x);

			size_t i;
			obj[0] = 0.;
			for (i = 0; i < m_num_vars; ++i) {
				obj[0] += fabs(x[i]);
			}

			obj[0] += m_bias;

			
			// evaluate constraint value

			for (auto &i : con)
				i = 0.;

			for (i = 0; i < m_num_vars; ++i)
			{
				con[0] += pow(x[i], 2);
			}
			con[0] -= 100 * m_num_vars;
			if (con[0] <= 0) con[0] = 0;

			con[1] = pow(cos(obj[0]) + sin(obj[0]), 2) - exp(cos(obj[0]) + sin(obj[0])) - 1 + exp(1.);
			if (fabs(con[1]) - 1e-4 <= 0) con[1] = 0;
			
		}
	}
}