#include "G14.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G14::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < m_num_vars; j++) {
				m_domain.setRange(0., 10., j);
			}
			for (size_t j = 0; j < m_num_vars; j++) {
				m_initial_domain.setRange(0., 10., j);
			}
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 3;
			m_constraint.resize(3);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kEquality;
			}

			//loadTranslation("/instance/problem/continuous/constrained/CEC2006/data/");  //data path 
			//setOriginalGlobalOpt();
			//setGlobalOpt(m_translation.data());
		}

		void G14::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			Real t = 0.0;
			Real c[10] = { -6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.1, -10.708, -26.662, -22.179 };
			for (int i = 0; i < 10; ++i)
			{
				t += x[i];
			}
			for (int i = 0; i < 10; ++i)
			{
				obj[0] += x[i] * (c[i] + log(x[i] /t));
			}

			// evaluate constraint value
			
			con[0] = x[0] + x[1] * 2 + x[2] * 2 + x[5] + x[9] - 2;
			con[1] = x[3] + 2 * x[4] + x[5] + x[6] - 1;
			con[2] = x[2] + x[6] + x[7] + 2 * x[8] + x[9] - 1;
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;


		}
	}
}

