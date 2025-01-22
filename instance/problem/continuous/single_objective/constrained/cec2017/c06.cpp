#include "C06.h"
#include "../../../../../core/instance_manager.h"

namespace ofec {
	namespace CEC2017 {
	
		void C06::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			setDomain(-20., 20.);
			setInitialDomain(-20., 20.);
			m_num_cons = 6;
			m_constraint.resize(6);
			for (auto& i : m_constraint)
				i = Constraint::kEquality;

			loadTranslation("/instance/problem/continuous/constrained/CEC2017/data/");  //data path
			setOriginalGlobalOpt(m_translation.data());
			m_optima = m_original_optima;
		}
		void C06::evaluateObjAndCon(Real *x, std::vector<Real>& obj, std::vector<Real> &con) {
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
			
			for (auto &i : con)
				i = 0.;

			for (i = 0; i < m_num_vars; ++i)
			{
				con[0] += x[i] * sin(x[i]);      //h1
			}
			con[0] *= -1;
			for (i = 0; i < m_num_vars; ++i)
			{
				con[1] += x[i] * sin(OFEC_PI*x[i]);   //h2
			}
			for (i = 0; i < m_num_vars; ++i)
			{
				con[2] += x[i] * cos(x[i]);    //h3
			}
			con[2] *= -1;
			for (i = 0; i < m_num_vars; ++i)
			{
				con[3] += x[i] * cos(OFEC_PI*x[i]);    //h4
			}
			for (i = 0; i < m_num_vars; ++i)
			{
				con[4] += x[i] * sin(2 * sqrt(fabs(x[i])));    //h5
			}
			for (i = 0; i < m_num_vars; ++i)
			{
				con[5] += x[i] * sin(2 * sqrt(fabs(x[i])));    //h6
			}
			con[5] *= -1;

			for (auto &i : con)
				if (fabs(i) - 1e-4 <= 0) i = 0;

		}
	}
}