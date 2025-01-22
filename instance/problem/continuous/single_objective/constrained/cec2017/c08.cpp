#include "C08.h"
#include "../../../../../core/instance_manager.h"
#include <math.h>

namespace ofec {
	namespace CEC2017 {
	
		void C08::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < m_num_vars; j++) {
				m_domain.setRange(-100., 100., j);
			}
			m_domain_update = true;
			for (size_t j = 0; j < m_num_vars; j++) {
				m_initial_domain.setRange(-100., 100., j);
			}
			m_num_cons = 2;
			m_constraint.resize(2);
			for (auto& i : m_constraint)
				i = Constraint::kEquality;

			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C08::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			obj[0] = *std::max_element(x, x + m_num_vars);

			obj[0] += m_bias;

			// evaluate constraint value
			
			std::vector<Real> y(x, x + m_num_vars);
			std::vector<Real> w(x, x + m_num_vars);

			for (auto &i : con) i = 0.;
			
			for (int a = 0; a < m_num_vars / 2; ++a) {
				y[a] = x[2 * a];
			}
			for (int b = 0; b < m_num_vars / 2; ++b)
			{
				Real m = 0.;
				for (int c = 0; c < b + 1; c++)    //h1
				{
					m += y[c];
				}
				con[0] += m*m;
			}

			for (int d = 0; d < m_num_vars / 2; ++d) {
				w[d] = x[2 * d+1];
			}
			for (int e = 0; e < m_num_vars / 2; ++e)
			{
				Real m = 0.;
				for (int f = 0; f < e + 1; f++)    //h2
				{
					m += w[f];
				}
				con[1] += m*m;
			}

			if (fabs(con[0]) - 1e-4 <= 0) con[0] = 0;
			if (fabs(con[1]) - 1e-4 <= 0) con[1] = 0;
			
		}
	}
}