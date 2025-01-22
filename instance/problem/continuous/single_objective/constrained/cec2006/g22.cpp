#include "G22.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G22::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(0., 20000., 0);
			for (size_t j = 1; j < 4; j++) {
				m_domain.setRange(0., 1000000., j);
			}
			for (size_t j = 4; j < 7; j++) {
				m_domain.setRange(0., 40000000., j);
			}
			m_domain.setRange(100., 299.99, 7);
			m_domain.setRange(100., 399.99, 8);
			m_domain.setRange(100.01, 300., 9);
			m_domain.setRange(100., 400., 10);
			m_domain.setRange(100., 600., 11);
			for (size_t j = 12; j < 15; j++) {
				m_domain.setRange(0., 500., j);
			}
			m_domain.setRange(0.01, 300., 15);
			m_domain.setRange(0.01, 400., 16);
			for (size_t j = 17; j < 22; j++) {
				m_domain.setRange(-4.7, 6.25, j);
			}

			m_initial_domain.setRange(0., 20000., 0);
			for (size_t j = 1; j < 4; j++) {
				m_initial_domain.setRange(0., 1000000., j);
			}
			for (size_t j = 4; j < 7; j++) {
				m_initial_domain.setRange(0., 40000000., j);
			}
			m_initial_domain.setRange(100., 299.99, 7);
			m_initial_domain.setRange(100., 399.99, 8);
			m_initial_domain.setRange(100.01, 300., 9);
			m_initial_domain.setRange(100., 400., 10);
			m_initial_domain.setRange(100., 600., 11);
			for (size_t j = 12; j < 15; j++) {
				m_initial_domain.setRange(0., 500., j);
			}
			m_initial_domain.setRange(0.01, 300., 15);
			m_initial_domain.setRange(0.01, 400., 16);
			for (size_t j = 17; j < 22; j++) {
				m_initial_domain.setRange(-4.7, 6.25, j);
			}
			//setDomain(0., 10., 0, 20);
			//setInitialDomain(0., 10., 0, 20);
			m_num_cons = 20;
			m_constraint.resize(20);
			m_constraint[0] = Constraint::kInequality;
			for (size_t i = 1; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kEquality;
			}

		}

		void G22::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {
			
			obj[0] = x[0];

			// evaluate constraint value
		
			con[0] =-x[0]+pow(x[1],0.6)+ pow(x[2], 0.6)+ pow(x[3], 0.6);
			con[1] = x[4] - 100000 * x[7] + 10000000;
			con[2] = x[5] + 100000 * x[7] - 100000 * x[8];
			con[3] = x[6] + 100000 * x[8] - 50000000;
			con[4] = x[4] + 100000 * x[9] - 33000000;
			con[5] = x[5] + 100000 * x[10] - 44000000;
			con[6] = x[6] + 100000 * x[11] - 66000000;
			con[7] = x[4] - 120 * x[1] * x[12];
			con[8] = x[5] - 80 * x[2] * x[13];
			con[9] = x[6] - 40 * x[3] * x[14];
			con[10] = x[7] - x[10] + x[15];
			con[11] = x[8] - x[11] + x[16];
			con[12] = -x[17] + log(x[9] - 100);
			con[13] = -x[18] + log(x[7] + 300);
			con[14] = -x[19] + log(x[15]);
			con[15] = -x[20] + log(-x[8] + 400);
			con[16] = -x[21] + log(x[16]);
			con[17] = -x[7] - x[9] + x[12] * x[17] - x[12] * x[18] + 400;
			con[18] = x[7] - x[8] - x[10] + x[13] * x[19] - x[13] * x[20] + 400;
			con[19] = x[8] - x[11] - 4.60517 * x[14] + x[14] * x[21] + 100;

			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3] <= 0) con[3] = 0;
			if (con[4] <= 0) con[4] = 0;
			if (con[5] <= 0) con[5] = 0;
			if (con[6] <= 0) con[6] = 0;
			if (con[7] <= 0) con[7] = 0;
			if (con[8] <= 0) con[8] = 0;
			if (con[9] <= 0) con[9] = 0;
			if (con[10] <= 0) con[10] = 0;
			if (con[11] <= 0) con[11] = 0;
			if (con[12] <= 0) con[12] = 0;
			if (con[13] <= 0) con[13] = 0;
			if (con[14] <= 0) con[14] = 0;
			if (con[15] <= 0) con[15] = 0;
			if (con[16] <= 0) con[16] = 0;
			if (con[17] <= 0) con[17] = 0;
			if (con[18] <= 0) con[18] = 0;
			if (con[19] <= 0) con[19] = 0;


		}
	}
}

