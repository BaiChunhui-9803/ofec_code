#include "C04.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {

		void C04::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-10., 10.);
			setInitialDomain(-10., 10.);
			m_num_cons = 2;
			m_constraint.resize(2);
			for (auto& i : m_constraint)
				i = Constraint::kInequality;

			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}

		void C04::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			size_t i;
			obj[0] = 0.;
			for (i = 0; i<m_num_vars; ++i)
			{
				obj[0] += (x[i] * x[i] - 10.0*cos(2.0*OFEC_PI*x[i]) + 10.0);
			}
			obj[0] += m_bias;

			
			// evaluate constraint value

			for (auto& i : con)
				i = 0.;
			for (i = 0; i < m_num_vars; ++i)
			{
				con[0] += x[i] * sin(2 * x[i]);
			}
			for (i = 0; i < m_num_vars; ++i)
			{
				con[1] += x[i] * sin(x[i]);
			}
			for (auto& i : con)
				if (i <= 0) i = 0;
		}
	}
}