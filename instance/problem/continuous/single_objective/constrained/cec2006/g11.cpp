#include "G11.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G11::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < m_num_vars; j++) {
				m_domain.setRange(-1., 1., j);
			}
			for (size_t j = 0; j < m_num_vars; j++) {
				m_initial_domain.setRange(-1., 1., j);
			}
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 1;
			m_constraint.resize(1);
			m_constraint[0] = Constraint::kEquality;

			//loadTranslation("/instance/problem/continuous/constrained/CEC2006/data/");  //data path 
			//setOriginalGlobalOpt();
			//setGlobalOpt(m_translation.data());
		}

		void G11::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			size_t i, j;
			obj[0] = x[0]*x[0] + (x[1] - 1)* (x[1] - 1);

			// evaluate constraint value
			
			con[0] = x[1] - x[0]*x[0] ;
			
			if (con[0] <= 0) con[0] = 0;


		}
	}
}

