
#include "ShiftedSchwefel_F15.h"
namespace ofec {
	namespace CEC2013 {
		ShiftedSchwefel_F15::ShiftedSchwefel_F15(const ParameterMap &v) : 
			ShiftedSchwefel_F15((v.at("problem name")), (v.at("number of variables")), 1) \
		{
			
		}

		ShiftedSchwefel_F15::ShiftedSchwefel_F15(const std::string &name, size_t size_var, size_t size_obj) : problem(name, size_var, size_obj), \
			function_CEC2013(name, size_var, size_obj) \
		{
			
		}

		ShiftedSchwefel_F15::~ShiftedSchwefel_F15() {
			delete[] mp_Ovector;

			delete[] mp_anotherz;

		}

		void ShiftedSchwefel_F15::initialize() {
			m_variable_monitor = true;
			setDomain(-100, 100);
			setInitialDomain(-100, 100);
			ID = 15;
			mp_anotherz = new Real[m_number_variables];

			// Ovector = createShiftVector(dimension,minX,maxX);
			mp_Ovector = readOvector();

			//setOriginalGlobalOpt();	
			setGlobalOpt(mp_Ovector);
			m_initialized = true;
		}

		int ShiftedSchwefel_F15::evaluateObjective(Real *x, std::vector<Real> &obj) {
			int i;
			Real result = 0.0;



			for (i = 0; i < m_number_variables; i++)
			{
				mp_anotherz[i] = x[i] - mp_Ovector[i];
			}

			result = schwefel(mp_anotherz, m_number_variables);

			obj[0] = result;
			return kNormalEval;
		}
	}
}