#include "C10.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {
	
		void C10::initialize_() {
			Function::initialize_();
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			m_num_cons = 2;
			m_constraint.resize(2);
			for (auto& i : m_constraint)
				i = Constraint::kEquality;


			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C10::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
			for (size_t i = 0; i < m_num_vars; ++i)
				x[i] -= m_translation[i];

			size_t i,j;
			obj[0] = *std::max_element(x, x + m_num_vars);

			obj[0] += m_bias;

			
			// evaluate constraint value

			for (auto &i : con)
				i = 0.;

			for (i = 1; i < m_num_vars + 1; ++i)
			{
				Real m = 0.;

				for (j = 0; j<i; j++)
				{
					m += x[j];
				}
				con[0] += m*m;
			}

			for (i = 0; i < m_num_vars - 1; ++i)
			{
				con[1] += pow((x[i] - x[i + 1]), 2);
			}

			for (auto &i : con)
				if (fabs(i) - 1e-4 <= 0) i = 0;
		}
	}
}