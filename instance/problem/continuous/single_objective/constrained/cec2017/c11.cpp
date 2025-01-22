#include "C11.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {
	
		void C11::initialize_() {
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
		void C11::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			size_t i;
			obj[0] = 0.;
			for (i = 0; i < m_num_vars ; ++i)
			{
				obj[0] += x[i];
			}

			obj[0] += m_bias;

			
			// evaluate constraint value

			for (auto &i : con)
				i = 0.;

			for (i = 0; i < m_num_vars; ++i)
			{
				con[0] *= x[i];
			}
			if (con[0] <= 0) con[0] = 0;

			for (i = 0; i < m_num_vars - 1; ++i)
			{
				con[1] += (x[i] - x[i + 1])*(x[i] - x[i + 1]);
			}
			if (fabs(con[1]) - 1e-4 <= 0) con[1] = 0;
	

		}
	}
}