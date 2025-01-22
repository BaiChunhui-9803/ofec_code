#include "G10.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G10::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(100., 10000., 0);
			for (size_t j = 1; j < 3; j++) {
				m_domain.setRange(1000., 10000., j);
			}
			for (size_t j = 3; j < m_num_vars; j++) {
				m_domain.setRange(10., 1000., j);
			}
			m_initial_domain.setRange(100., 10000., 0);
			for (size_t j = 1; j < 3; j++) {
				m_initial_domain.setRange(1000., 10000., j);
			}
			for (size_t j = 3; j < m_num_vars; j++) {
				m_initial_domain.setRange(10., 1000., j);
			}
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 2;
			m_constraint.resize(2);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}

			//loadTranslation("/instance/problem/continuous/constrained/CEC2006/data/");  //data path 
			//setOriginalGlobalOpt();
			//setGlobalOpt(m_translation.data());
		}

		void G10::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = x[0] + x[1] + x[2];

			// evaluate constraint value
			
			con[0] = -1 + 0.0025 * (x[3] + x[5]);
			con[1] = -1 + 0.0025 * (x[4] + x[6] - x[3]);
			con[2] = -1 + 0.01 * (x[7] - x[4]);
			con[3] = -x[0] * x[5] + 833.33252 * x[3] + 100 * x[0] - 83333.333;
			con[4] = -x[1] * x[6] + 1250 * x[4] + x[1] * x[3] - 1250 * x[3];
			con[5] = -x[2] * x[7] + 1250000 + x[2] * x[4] - 2500 * x[4];
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;
			if (con[5] <= 0) con[5] = 0;


		}
	}
}

