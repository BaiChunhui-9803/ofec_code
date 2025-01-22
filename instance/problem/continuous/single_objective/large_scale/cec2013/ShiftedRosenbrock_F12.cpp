
#include "ShiftedRosenbrock_F12.h"
namespace ofec {
	namespace CEC2013 {
		ShiftedRosenbrock_F12::ShiftedRosenbrock_F12(const ParameterMap &v) : 
			ShiftedRosenbrock_F12((v.at("problem name")), (v.at("number of variables")), 1) \
		{
			
		}

		ShiftedRosenbrock_F12::ShiftedRosenbrock_F12(const std::string &name, size_t size_var, size_t size_obj) : problem(name, size_var, size_obj), \
			function_CEC2013(name, size_var, size_obj) \
		{
			
		}

		ShiftedRosenbrock_F12::~ShiftedRosenbrock_F12() {
			delete[] mp_Ovector;

			delete[] mp_anotherz;

		}

		void ShiftedRosenbrock_F12::initialize() {
			m_variable_monitor = true;
			setDomain(-100, 100);
			setInitialDomain(-100, 100);
			ID = 12;
			mp_anotherz = new Real[m_number_variables];

			// Ovector = createShiftVector(dimension,minX,maxX);
			mp_Ovector = readOvector();

			//setOriginalGlobalOpt();
			std::vector<Real> data(m_number_variables,1);
			for (size_t i = 0; i < m_number_variables; ++i)
				data[i] += mp_Ovector[i];
			setGlobalOpt(data.data());
			m_initialized = true;
		}

		int ShiftedRosenbrock_F12::evaluateObjective(Real *x, std::vector<Real> &obj) {
			int i;
			Real result = 0.0;

			for (i = 0; i < m_number_variables; i++)
			{
				mp_anotherz[i] = x[i] - mp_Ovector[i];
			}

			result = rosenbrock(mp_anotherz, m_number_variables);

			obj[0] = result;
			return kNormalEval;
		}
	}
}