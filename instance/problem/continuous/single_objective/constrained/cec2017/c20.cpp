#include "C20.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {

		void C20::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			m_num_cons = 2;
			m_constraint.resize(2);
			for (auto& i : m_constraint)
				i = Constraint::kInequality;

			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C20::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];
			

			size_t i;
			obj[0] = 0.;

			for (i = 0; i<m_num_vars; ++i)
			{
				obj[0] += 0.5 + (pow(sin(pow((x[i] * x[i] + x[i + 1] * x[i + 1]), 0.5)), 2) - 0.5) / pow(1 + 0.001*pow(x[i] * x[i] + x[i + 1] * x[i + 1], 0.5), 2);
			}
			auto n = m_num_vars-1;
			obj[0] += 0.5 + (pow(sin(pow((x[n] * x[n] + x[0] * x[0]), 0.5)), 2) - 0.5) / pow(1 + 0.001*pow(x[n] * x[n] + x[0] * x[0], 0.5), 2);
			obj[0] += m_bias;

			// evaluate constraint value
			
			for (auto &i : con)
				i = 0.;

			Real num1 = 0.;
			for (i = 0; i < m_num_vars; ++i)
				num1 += x[i];

			con[0] = pow(cos(num1), 2) - 0.25*cos(num1) - 0.125;
			if (con[0] <= 0) con[0] = 0;

			con[1] = exp(cos(num1)) - exp(0.25);
			if (con[1] <= 0) con[1] = 0;

		}
	}
}