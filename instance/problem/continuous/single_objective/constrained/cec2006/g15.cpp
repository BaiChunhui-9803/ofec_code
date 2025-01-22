#include "G15.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G15::initialize_() {
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
			m_num_cons = 2;
			m_constraint.resize(2);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kEquality;
			}

			//loadTranslation("/instance/problem/continuous/constrained/CEC2006/data/");  //data path 
			//setOriginalGlobalOpt();
			//setGlobalOpt(m_translation.data());
		}

		void G15::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {
			obj[0] = 1000 - x[0] * x[0] - 2 * x[1] * x[1] - x[2] * x[2] - x[0] * x[1] - x[0] * x[2];

			// evaluate constraint value
			Real t1 = 1.0, t2 = 0.0;
			con[0] = x[0] *x[0] + x[1] * x[1] + x[2] * x[2] - 25;
			con[1] = 8 * x[0] + 14 * x[1] + 7 * x[2] - 56;
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;


		}
	}
}

