#include "G05.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G05::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < 2; j++) {
				m_domain.setRange(0., 1200., j);
			}
			for (size_t j = 2; j < 4; j++) {
				m_domain.setRange(-0.55, 0.55, j);
			}
			m_domain_update = true;
			for (size_t j = 0; j < 2; j++) {
				m_initial_domain.setRange(0., 1200., j);
			}
			for (size_t j = 2; j < 4; j++) {
				m_initial_domain.setRange(-0.55, 0.55, j);
			}
			m_num_cons = 5;
			m_constraint.resize(5);
			for (size_t i = 0; i < 2; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}
			for (size_t i = 2; i < 5; ++i)
			{
				m_constraint[i] = Constraint::kEquality;
			}

			//loadTranslation("/instance/problem/continuous/constrained/CEC2006/data/");  //data path 
			//setOriginalGlobalOpt();
			//setGlobalOpt(m_translation.data());
		}

		void G05::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = 3 * x[0] + 0.000001 * x[0] * x[0] * x[0] + 2 * x[1] + (0.000002 / 3) * x[1] * x[1] * x[1];

			// evaluate constraint value
			con[0] = -x[3] + x[2] - 0.55;
			con[1] = -x[2] + x[3] - 0.55;
			con[2] = 1000 * sin(-x[2] - 0.25) + 1000 *sin(-x[3] - 0.25) + 894.8 - x[0];
			con[3] = 1000 * sin(x[2] - 0.25) + 1000 * sin(x[2] - x[3] - 0.25) + 894.8 - x[1];
			con[4] = 1000 * sin(x[3] - 0.25) + 1000 * sin(x[3] - x[2] - 0.25) + 1294.8;
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;


		}
	}
}

