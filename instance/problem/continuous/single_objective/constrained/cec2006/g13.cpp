#include "G13.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G13::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < 2; j++) {
				m_domain.setRange(-2.3, 2.3, j);
			}
			for (size_t j = 2; j < m_num_vars; j++) {
				m_domain.setRange(-3.2, 3.2, j);
			}
			for (size_t j = 0; j < 2; j++) {
				m_initial_domain.setRange(-2.3, 2.3, j);
			}
			for (size_t j = 2; j < m_num_vars; j++) {
				m_initial_domain.setRange(-3.2, 3.2, j);
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

		void G13::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = exp(x[0] * x[1] * x[2] * x[3] * x[4]);

			// evaluate constraint value
			con[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4] - 10;
			con[1] = x[1] * x[2] - 5 * x[3] * x[4];
			con[2] = x[0] * x[0] + x[1] * x[1] + 1;
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;


		}
	}
}

