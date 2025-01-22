#include "G16.h"
#include "../../../../../core/instance_manager.h"
#include <cmath>


namespace ofec {
	namespace CEC2006 {

		void G16::initialize_() {
			Function::initialize_();
			//m_variable_monitor = true;
			auto& v = GET_PARAM(m_id_param);
			resizeVariable(v.get<int>("number of variables"));
			m_domain.setRange(704.4148, 906.3855, 0);
			m_domain.setRange(68.6, 288.88, 1);
			m_domain.setRange(0., 134.75, 2);
			m_domain.setRange(193, 287.0966, 3);
			m_domain.setRange(25, 84.1988, 4);
			m_initial_domain.setRange(704.4148, 906.3855, 0);
			m_initial_domain.setRange(68.6, 288.88, 1);
			m_initial_domain.setRange(0., 134.75, 2);
			m_initial_domain.setRange(193, 287.0966, 3);
			m_initial_domain.setRange(25, 84.1988, 4);

			m_num_cons = 38;
			m_constraint.resize(38);
			for (size_t i = 0; i < m_num_cons; ++i)
			{
				m_constraint[i] = Constraint::kInequality;
			}

			//loadTranslation("/instance/problem/continuous/constrained/CEC2006/data/");  //data path 
			//setOriginalGlobalOpt();
			//setGlobalOpt(m_translation.data());
		}

		void G16::evaluateObjAndCon(Real* x, std::vector<Real>& obj, std::vector<Real>& con) {

			Real y1 = x[1] + x[2] + 41.6;
			Real c1 = 0.024 * x[3] - 4.62;
			Real y2 = 12.5 / c1 + 12;
			Real c2 = 0.0003535 * x[0] * x[0] + 0.5311 * x[0] + 0.08705 * y2 * x[0];
			Real c3 = 0.052 * x[0] + 78 + 0.002377*y2 * x[0];
			Real y3 = c2 / c3;
			Real y4 = 19 * y3;
			Real c4 = 0.04782 * (x[0] - y3) + 0.1956 * (x[0] - y2) * (x[0] - y2) / x[1] + 0.6376 * y3 + 1.594 * y3;
			Real c5 = 100 * x[1];
			Real c6 = x[0] - y3 - y4;
			Real c7 = 0.95 - c4 / c5;
			Real y5 = c6 * c7;
			Real y6 = x[0] - y5 - y4 - y3;
			Real c8 = (y5 + y4) * 0.995;
			Real y7 = c8 / y1;
			Real y8 = c8 / 3798;
			Real c9 = y7 - 0.0663 * y7 / y8 - 0.3135;
			Real y9 = 96.82 / c9 + 0.321 * y1;
			Real y10 = 1.29 * y5 + 1.258 * y4 + 2.29 * y3 + 1.71 * y6;
			Real y11 = 1.71 * x[0] - 0.452 * y4 + 0.58 * y3;
			Real c10 = 12.3 / 752.3;
			Real c11 = (1.75 * y2) * (0.995 * x[0]);
			Real c12 = 0.995 * y10 + 1998;
			Real y12 = c10 * x[0] + c11 / c12;
			Real y13 = c12 - 1.75 * y2;
			Real y14 = 3623 + 64.4 * x[1] + 58.4 * x[2] + 146312 / (y9 + x[4]);
			Real c13 = 0.995 * y10 + 60.8 * x[1] + 48 * x[3] - 0.1121 * y14 - 5095;
			Real y15 = y13 / c13;
			Real y16 = 148000 - 331000 * y15 + 40 * y13 - 61 * y15 * y13;
			Real c14 = 2324 * y10 - 28740000 * y2;
			Real y17 = 14130000 - 1328 * y10 - 531 * y11 + c14 / c12;
			Real c15 = y13 / y15 - y13 / 0.52;
			Real c16 = 1.104 - 0.72 * y15;
			Real c17 = y9 + x[4];

			obj[0] = 0.000117 * y14 + 0.1365 + 0.00002358 * y13 + 0.000001502 * y16 + 0.0321 * y12 + 0.004324 * y5 + 0.0001 * c15 / c16 + 37.48 * y2 / c12 - 0.0000005843 * y17;

			// evaluate constraint value
			con[0] = 0.28 / 0.72 * y5 - y4;
			con[1] = x[2] - 1.5 * x[1];
			con[2] = 3496 * y2 / c12 - 21;
			con[3] = 110.6 + y1 - 62212 / c17;
			con[4] = 213.1 - y1;
			con[5] = y1 - 405.23;
			con[6] = 17.505 - y2;
			con[7] = y2 - 1053.6667;
			con[8] = 11.275 - y3;
			con[9] = y3 - 35.03;
			con[10] = 214.228 - y4;
			con[11] = y4 - 665.585;
			con[12] = 7.458 - y5;
			con[13] = y5 - 584.463;
			con[14] = 0.961 - y6;
			con[15] = y6 - 265.916;
			con[16] = 1.612 - y7;
			con[17] = y7 - 7.046;
			con[18] = 0.146 - y8;
			con[19] = y8 - 0.222;
			con[20] = 107.99 - y9;
			con[21] = y9 - 273.366;
			con[22] = 922.693 - y10;
			con[23] = y10 - 1286.105;
			con[24] = 926.832 - y11;
			con[25] = y11 - 1444.046;
			con[26] = 18.766 - y12;
			con[27] = y12 - 537.141;
			con[28] = 1072.163 - y13;
			con[29] = y13 - 3247.039;
			con[30] = 8961.448 - y14;
			con[31] = y14 - 26844.086;
			con[32] = 0.063 - y15;
			con[33] = y15 - 0.386;
			con[34] = 71084.33 - y16;
			con[35] = -140000 + y16;
			con[36] = 2802713 - y17;
			con[37] = y17 - 12146108;

			if (con[0] <= 0) con[0] = 0;
			if (con[1] <= 0) con[1] = 0;
			if (con[2] <= 0) con[2] = 0;
			if (con[3]<= 0) con[3] = 0;
			if (con[4]<= 0) con[4] = 0;
			if (con[5]<= 0) con[5] = 0;
			if (con[6] <= 0) con[6] = 0;
			if (con[7] <= 0) con[7] = 0;
			if (con[8] <= 0) con[8] = 0;
			if (con[9] <= 0) con[9] = 0;
			if (con[10]<= 0) con[10]= 0;
			if (con[11]<= 0) con[11]= 0;
			if (con[12]<= 0) con[12]= 0;
			if (con[13]<= 0) con[13]= 0;
			if (con[14]<= 0) con[14]= 0;
			if (con[15]<= 0) con[15]= 0;
			if (con[16]<= 0) con[16]= 0;
			if (con[17]<= 0) con[17]= 0;
			if (con[18]<= 0) con[18]= 0;
			if (con[19]<= 0) con[19]= 0;
			if (con[20]<= 0) con[20]= 0;
			if (con[21]<= 0) con[21]= 0;
			if (con[22]<= 0) con[22]= 0;
			if (con[23]<= 0) con[23]= 0;
			if (con[24]<= 0) con[24]= 0;
			if (con[25]<= 0) con[25]= 0;
			if (con[26]<= 0) con[26]= 0;
			if (con[27]<= 0) con[27]= 0;
			if (con[28] <= 0) con[28]= 0;
			if (con[29] <= 0) con[29]= 0;
			if (con[30] <= 0) con[30]= 0;
			if (con[31] <= 0) con[31]= 0;
			if (con[32] <= 0) con[32]= 0;
			if (con[33] <= 0) con[33]= 0;
			if (con[34] <= 0) con[34]= 0;
			if (con[35] <= 0) con[35]= 0;
			if (con[36] <= 0) con[36]= 0;
			if (con[37] <= 0) con[37]= 0;


		}
	}
}

