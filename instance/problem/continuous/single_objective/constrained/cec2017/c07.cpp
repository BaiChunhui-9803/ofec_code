#include "C07.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {
		
		void C07::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-50., 50.);
			setInitialDomain(-50., 50.);
			m_num_cons = 2;
			m_constraint.resize(2);
			for (auto& i : m_constraint)
				i = Constraint::kEquality;

			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C07::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			size_t i;
			obj[0] = 0.;
			for (i = 0; i<m_num_vars - 1; ++i)
			{
				obj[0] += x[i] * sin(x[i]);
			}
			obj[0] += m_bias;

			
			// evaluate constraint value
			
			for (auto &i : con)
				i = 0.;

			for (i = 0; i < m_num_vars; ++i)
			{
				con[0] += (x[i] - 100 * cos(0.5*x[i]) + 100);       
			}
			
			con[1] *= -1 * con[0];

			for (auto &i : con)
				if (fabs(i) - 1e-4 <= 0) i = 0;
		}
	}
}