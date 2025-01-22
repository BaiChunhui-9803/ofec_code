#include "C09.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {
	
		void C09::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-10., 10.);
			setInitialDomain(-10., 10.);
			m_num_cons = 2;
			m_constraint.resize(2);
			m_constraint[0] = Constraint::kInequality;
			m_constraint[1] = Constraint::kEquality;


			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C09::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			size_t i;
			obj[0] = *std::max_element(x, x + m_num_vars);

			obj[0] += m_bias;

			// evaluate constraint value

			std::vector<Real> x_1(x, x + m_num_vars);
			std::vector<Real> x_2(x, x + m_num_vars);

			for (auto &i : con)
				i = 0.;

			for (i = 1; i <= m_num_vars / 2; ++i) {
				x_2[i - 1] = x[2 * i - 1];
			}
			for (i = 0; i < m_num_vars / 2; ++i)
			{
				con[0] *= x_2[i];
			}
			if (con[0] <= 0) con[0] = 0;

			for (i = 1; i <= m_num_vars / 2; ++i) {
				x_1[i - 1] = x[2 * i - 2];
			}
			for (i = 0; i < m_num_vars / 2; ++i)
			{
				con[1] += (x_1[i] * x_1[i] - x_1[i + 1])*(x_1[i] * x_1[i] - x_1[i + 1]);
			}
			if (fabs(con[1]) - 1e-4 <= 0) con[1] = 0;
		}
	}
}