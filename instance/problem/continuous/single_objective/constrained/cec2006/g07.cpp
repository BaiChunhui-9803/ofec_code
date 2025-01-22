#include "G07.h"
#include "../../../../../core/instance_manager.h"


namespace ofec {
	namespace CEC2006 {

		void G07::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < m_num_vars; j++) {
				m_domain.setRange(-10., 10., j);
			}
			m_domain_update = true;
			for (size_t j = 0; j < m_num_vars; j++) {
				m_initial_domain.setRange(-10., 10., j);
			}
			m_num_cons = 8;
			m_constraint.resize(8);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}

			//loadTranslation("/instance/problem/continuous/constrained/CEC2006/data/");  //data path 
			//setOriginalGlobalOpt();
			//setGlobalOpt(m_translation.data());
		}

		void G07::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = x[0] * x[0] + x[1] * x[1] + x[0] * x[1] - 14 * x[0] - 16 * x[1] + (x[2] - 10) * (x[2] - 10) + 4 * (x[3] - 5) * (x[3] - 5) + (x[4] - 3) * (x[4] - 3) + 2 * (x[5] - 1) * (x[5] - 1) + 5 * x[6] * x[6] + 7 * (x[7] - 11) * (x[7] - 11) + 2 * (x[8] - 10) * (x[8] - 10) + (x[9] - 7) * (x[9] - 7) + 45;

			// evaluate constraint value

			con[0] = -105 + 4 * x[0] + 5 * x[1] - 3 * x[6] + 9 * x[7];
			con[1] = 10 * x[0] - 8 * x[1] - 17 * x[6] + 2 * x[7];
			con[2] = -8 * x[0] + 2 * x[1] + 5 * x[8] - 2 * x[9] - 12;
			con[3] = 3 * (x[0] - 2) * (x[0] - 2) + 4 * (x[1] - 3) * (x[1] - 3) + 2 * x[2] * x[2] - 7 * x[3] - 120;
			con[4] = 5 * x[0] * x[0] + 8 * x[1] + (x[2] - 6) * (x[2] - 6) - 2 * x[3] - 40;
			con[5] = x[0] * x[0] + 2 * (x[1] - 2) * (x[1] - 2) - 2 * x[0] * x[1] + 14 * x[4] - 6 * x[5];
			con[6] = 0.5 * (x[0] - 8) * (x[0] - 8) + 2 * (x[1] - 4) * (x[1] - 4) + 3 * x[4] * x[4] - x[5] - 30;
			con[7] = -3 * x[0] + 6 * x[1] + 12 * (x[8] - 8) * (x[8] - 8) - 7 * x[9];
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;
			if (con[5] <= 0) con[5] = 0;
			if (con[6] <= 0) con[6] = 0;
			if (con[7] <= 0) con[7] = 0;

		}
	}
}

