#include "ATF1.h"

#include <math.h>


namespace ofec {
	namespace ATF {

		void _1::addInputParameters() {
			m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
		}
		void _1::initialize_(Environment* env) {
			Function::initialize_(env);
			//m_variable_monitor = true;
			//auto& v = GET_PARAM(m_id_param);
			//resizeVariable(v.get<int>("number of variables"));
			for (size_t j = 0; j < m_number_variables; j++) {
				m_domain.setRange(-5., 5., j);
			}
			//m_domain_update = true;
			//for (size_t j = 0; j < m_number_variables; j++) {
			//	m_initial_domain.setRange(-5., 5., j);
			//}
			m_number_constraints = 3;
			m_constraint.resize(3);
            for(int i=0;i<3;i++){
                m_constraint[i] = Constraint::kInequality;
            }
		}

		void _1::evaluateObjectiveAndConstraint(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {
			
			obj[0] = x[0]*x[0]-10*cos(2*3.14*x[0])+x[1]*x[1]-10*cos(2*3.14*x[1])+20;

			con[0] = x[0]+2;
            con[1] = -x[0]-3;
            con[2] = 2*x[0]+x[1]+7;
			
			if (con[0] <= 0) con[0] = 0;
            if (con[1] <= 0) con[1] = 0;
            if (fabs(con[2]) - 1e-4 <= 0) con[2] = 0;
		
		}
	}
}