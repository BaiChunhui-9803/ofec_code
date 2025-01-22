#include "cec2015_function.h"



#include "../../../../../../core/global.h"
#include "../../../../../../core/problem/solution.h"

namespace ofec {
	namespace CEC2015 {
		//CEC2015_function::CEC2015_function(const std::string &name, size_t size_var, size_t size_obj) : problem(name, size_var, size_obj), 
		//	function(name, size_var, size_obj)
		//{
		//	setDomain(-100, 100);
		//	setInitialDomain(-100, 100);
		//	
		//	m_variable_accuracy = 0.01;
		//	m_objective_accuracy = 1.e-4;
		//	setConditionNumber(1.0);
		//	m_variable_monitor = true;
		//}

		void CEC2015_function::initialize_()
		{
			Function::initialize_(env);
			resizeObjective(1);
			m_optimize_mode[0] = OptimizeMode::kMinimize;
			auto& v = *m_param;

			resizeVariable(v.get<int>("number of variables"));	
			//resizeVariable(v.get<int>("number of variables"));


			setDomain(-100, 100);
			setInitialDomain(-100, 100);

			m_variable_accuracy = 0.01;
			m_objective_accuracy = 1.e-4;
			setConditionNumber(1.0);
		//	m_variable_monitor = true;
		}
		
		bool CEC2015_function::loadTranslation(const std::string &path) {
			std::string s;
			std::stringstream ss;
			ss << ".txt";
			s = ss.str();
			s.insert(0, m_name + "_Shift");
			s.insert(0, path);    // data path
			s.insert(0, g_working_directory);

			load_translation_(s);

			return true;
		}
		void CEC2015_function::load_translation_(const std::string &path) {
			m_translation.resize(m_number_variables);
			std::ifstream in(path);
			if (in.fail()) {
				
				std::vector<Real> temp(m_number_variables, 20);
				setTranslation(temp.data());
				
				std::ofstream out(path);
				for (auto &i : m_translation)
					out << i << " ";
				out.close();
			}
			else {
				for (auto &i : m_translation)
					in >> i;
			}
			in.close();
			m_translated = true;
		}

		
		bool CEC2015_function::load_optima(const std::string &path) {
			std::string s;
			std::stringstream ss;
			ss << m_number_variables << "Dim.txt";
			s = ss.str();
			s.insert(0, m_name + "_Optima_");

			s.insert(0, path);// data path
			s.insert(0, g_working_directory);

			load_optima_(s);

			return true;
		}

		void CEC2015_function::load_optima_(const std::string &path) {
			std::ifstream in;
			in.open(path);
			if (in.fail()) {
				set_optima();
				std::ofstream out(path);
				for (size_t i = 0; i < m_optima->numberSolutions();++i)
					for (size_t j = 0; j < m_number_variables; ++j)
						out << optima()->solution(i).variable()[j] << " ";
				out.close();
			}
			else {
				while(!(in.eof())) {
					VariableVector<Real> temp_var(m_number_variables);
					for (size_t j = 0; j < m_number_variables; ++j)
						in >> temp_var[j];
					optima()->appendVar(std::move(temp_var));
				}
			}
			in.close();
		}

		void CEC2015_function::set_optima() {
			//  to do ..
			//throw MyExcept("Waiting to do !");
		}

		void CEC2015_function::evaluate_optima(){
            Solution<VariableVector<Real>> s(m_number_objectives, m_number_constraints, m_number_variables);
            for (size_t i = 0; i < m_optima->numberSolutions(); ++i) {
			    s.variable() = optima()->solution(i).variable();
				evaluate_(s,false);
				//s.evaluate(this, -1, false);
			   // s.evaluate(false, caller::Problem);
				m_optima->appendObj(s.objective());
			}
			
		}
		void CEC2015_function::rotate(Real *x) {
			Real *x_ = new Real[m_number_variables];
			std::copy(x, x + m_number_variables, x_);

			for (size_t i = 0; i<m_number_variables; ++i) {
				x[i] = 0;

				for (size_t j = 0; j < m_number_variables; ++j) {
					x[i] += m_rotation[i][j] * x_[j];
				}
			}

			delete[] x_;
			x_ = 0;
		}
	}
}