#include "N20_ShiftedSchwefel_F11.h"

namespace ofec {
	namespace CEC2013 {
		N20_ShiftedSchwefel_F11::N20_ShiftedSchwefel_F11(const ParameterMap &v) : 
			N20_ShiftedSchwefel_F11((v.at("problem name")), (v.at("number of variables")), 1) \
		{
			
		}

		N20_ShiftedSchwefel_F11::N20_ShiftedSchwefel_F11(const std::string &name, size_t size_var, size_t size_obj) : problem(name, size_var, size_obj), \
			function_CEC2013(name, size_var, size_obj) \
		{
			
		}

		N20_ShiftedSchwefel_F11::~N20_ShiftedSchwefel_F11() {
			delete[] mp_Ovector;
			delete[] mp_Pvector;
			delete[] mp_anotherz;

			for (int i = 0; i < 25; ++i)
				delete[] mpp_r25[i];
			for (int i = 0; i < 50; ++i)
				delete[] mpp_r50[i];
			for (int i = 0; i < 100; ++i)
				delete[] mpp_r100[i];

			delete[] mpp_r25;
			delete[] mpp_r50;
			delete[] mpp_r100;
			delete[] mp_s;
			delete[] mp_w;

		}

		void N20_ShiftedSchwefel_F11::initialize() {
			m_variable_monitor = true;
			setDomain(-32, 32);
			setInitialDomain(-32, 32);
			ID = 11;
			m_nonSeparableGroupNumber = 20;
			mp_anotherz = new Real[m_number_variables];

			// Ovector = createShiftVector(dimension,minX,maxX);
			mp_Ovector = readOvector();
			mp_Pvector = readPermVector();
			mpp_r25 = readR(25);
			mpp_r50 = readR(50);
			mpp_r100 = readR(100);
			mp_s = readS(m_nonSeparableGroupNumber);
			mp_w = readW(m_nonSeparableGroupNumber);

			//setOriginalGlobalOpt();
			setGlobalOpt(mp_Ovector);
			m_initialized = true;
		}

		int N20_ShiftedSchwefel_F11::evaluateObjective(Real *x, std::vector<Real> &obj) {
			size_t i;
			Real result = 0.0;

			for (i = 0; i < m_number_variables; i++)
			{
				mp_anotherz[i] = x[i] - mp_Ovector[i];
			}

			// s_size non-separable part with rotation
			size_t c = 0;
			for (i = 0; i < m_nonSeparableGroupNumber; i++)
			{
				mp_anotherz1 = rotate_vector(i, c);
				result += mp_w[i] * schwefel(mp_anotherz1, mp_s[i]);
				delete[] mp_anotherz1;
			}

			obj[0] = result;
			return kNormalEval;
		}
	}
}