#include "G04.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G04::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(78., 102., 0);
			m_domain.setRange(33., 45., 1);
			for (size_t j = 2; j < m_num_vars; j++) {
				m_domain.setRange(27., 45., j);
			}
			m_domain_update = true;
			m_initial_domain.setRange(78., 102., 0);
			m_initial_domain.setRange(33., 45., 1);
			for (size_t j = 2; j < m_num_vars; j++) {
				m_initial_domain.setRange(27., 45., j);
			}
			m_num_cons = 6;
			m_constraint.resize(6);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}

			//loadTranslation("/instance/problem/continuous/constrained/CEC2006/data/");  //data path 
			//setOriginalGlobalOpt();
			//setGlobalOpt(m_translation.data());
		}

		void G04::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = 5.3578547*x[2]*x[2]+0.8356891*x[0]*x[4]+37.293239*x[0]-40792.141;

			// evaluate constraint value
			
			con[0] = 85.334407 + 0.0056858 * x[1] * x[4] + 0.0006262 * x[0] * x[3] - 0.0022053 * x[2] * x[4] - 92;
			con[1] = -85.334407 - 0.0056858 * x[1] * x[4] - 0.0006262 * x[0] * x[3] + 0.0022053 * x[2] * x[4];
			con[2] = 80.51249 + 0.0071317 * x[1] * x[4] + 0.0029955 * x[0] * x[1] + 0.0021813 * x[2] * x[2] - 110;
			con[3] = -80.51249 - 0.0071317 * x[1] * x[4] - 0.0029955 * x[0] * x[1] - 0.0021813 * x[2] * x[2] + 90;
			con[4] = 9.300961 + 0.0047026 * x[2] * x[4] + 0.0012547 * x[0] * x[2] + 0.0019085 * x[2] * x[3] - 25;
			con[5] = -9.300961 - 0.0047026 * x[2] * x[4] - 0.0012547 * x[0] * x[2] - 0.0019085 * x[2] * x[3] + 20;
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;

		}
	}
}

