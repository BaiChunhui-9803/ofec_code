#include "G09.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G09::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < m_num_vars; j++) {
				m_domain.setRange(-10., 10., j);
			}
			for (size_t j = 0; j < m_num_vars; j++) {
				m_initial_domain.setRange(-10., 10., j);
			}
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 4;
			m_constraint.resize(4);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}

			//loadTranslation("/instance/problem/continuous/constrained/CEC2006/data/");  //data path 
			//setOriginalGlobalOpt();
			//setGlobalOpt(m_translation.data());
		}

		void G09::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			obj[0] = pow((x[0] - 10),2) + 5 * pow((x[1] - 12) ,2) + pow(x[2] ,4) + 3 * pow((x[3] - 11) ,2) + 10 * pow(x[4] ,6) + 7 * pow(x[5] ,2) + pow(x[6] ,4) - 4 * x[5] * x[6] - 10 * x[5] - 8 * x[6];

			// evaluate constraint value
			
			con[0] = -127 + 2 * (pow(x[0] ,2) + 3 * pow(x[1] ,4) + x[2] + 4 * pow(x[3] ,2) + 5 * x[4]);
			con[1] = -282 + 7 * x[0] + 3 * x[1] + 10 * pow(x[2] ,2) + x[3] - x[4];
			con[2] = -196 + 23 * x[0] + pow(x[1] ,2 )+ 6 * pow(x[5] ,2) - 8 * x[6];
			con[3] = 4 * pow(x[0] ,2) + pow(x[1] ,2) - 3 * x[0] * x[1] + 2 * pow(x[2] ,2) + 5 * x[5] - 11 * x[6];
			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;

		}
	}
}

