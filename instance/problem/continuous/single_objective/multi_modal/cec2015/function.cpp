#include "function.h"
#include <sstream>
#include <fstream>
#include "../../../../../../core/global.h"


namespace ofec::cec2015 {
	void Function::initialize_(Environment *env) {
		Continuous::initialize_();
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		int num_vars = m_param->get<int>("number of variables");
		resizeVariable(num_vars);
		setDomain(-100, 100);
		m_variable_niche_radius = 0.01;
		m_objective_accuracy = 1.e-4;
	}

	bool Function::loadOptima(const std::string &path) {
		std::string s;
		std::stringstream ss;
		ss << m_number_variables << "Dim.txt";
		s = ss.str();
		s.insert(0, m_name + "_Optima_");
		s.insert(0, path);// data path
		s.insert(0, g_working_directory);
		m_optima.reset(new Optima<>());
		loadOptima_(s);
		m_optima->setVariableGiven(true);
		return true;
	}

	void Function::loadOptima_(const std::string &path) {
		std::ifstream in;
		in.open(path);
		if (!in.fail()) {
			while (!(in.eof())) {
				VariableVector<Real> temp_var(m_number_variables);
				for (size_t j = 0; j < m_number_variables; ++j)
					in >> temp_var[j];
				optima()->appendVar(std::move(temp_var));
			}
		}
		else
		//	throw MyExcept("No file for optima data exists@cec2015::Function.");
		in.close();
	}


	bool Function::loadTranslation(const std::string &path) {
		std::string s;
		std::stringstream ss;
		ss << ".txt";
		s = ss.str();
		s.insert(0, m_name + "_Shift");
		s.insert(0, path);    // data path
		s.insert(0, g_working_directory);
		loadTranslation_(s);
		return true;
	}

	void Function::loadTranslation_(const std::string &path) {
		m_translation.resize(m_number_variables);
		std::ifstream in(path);
		if (!in.fail()) {
			for (auto &i : m_translation)
				in >> i;
		}
		else
			throw MyExcept("No file for tranlation data exists@cec2015::Function.");
		in.close();
	}

	bool Function::loadRotation(const std::string &path) {
		std::string s;
		std::stringstream ss;
		ss << m_number_variables << "Dim.txt";
		s = ss.str();
		s.insert(0, m_name + "_RotM_");
		s.insert(0, path);    // data path
		s.insert(0, g_working_directory);
		loadRotation_(s);
		return true;
	}

	void Function::loadRotation_(const std::string &path) {
		m_rotation.resize(m_number_variables);
		for (auto &row : m_rotation)
			row.resize(m_number_variables);
		std::ifstream in(path);
		if (!in.fail()) {
			for (auto &row : m_rotation) {
				for (size_t i = 0; i < m_number_variables; ++i)
					in >> row[i];
			}
		}
		else
			throw MyExcept("No file for rotation data exists@cec2015::Function.");
		in.close();
	}

	void Function::evaluateOptima() {
		if (m_optima->numberSolutions() > 0) {
			std::vector<Real> x(m_number_variables), obj(1);
			for (size_t i = 0; i < m_optima->numberSolutions(); ++i) {
				x = optima()->solution(i).variable().vect();
				evaluateObjective(x.data(), obj);
				m_optima->appendObj(obj);
			}
			m_optima->setObjectiveGiven(true);
		}
	}

	void Function::translate(Real *x) {
		for (size_t i = 0; i < m_number_variables; ++i)
			x[i] -= m_translation[i];
	}

	void Function::scale(Real *x) {
		for (size_t i = 0; i < m_number_variables; ++i)
			x[i] /= m_scale;
	}

	void Function::rotate(Real *x) {
		Real *x_ = new Real[m_number_variables];
		std::copy(x, x + m_number_variables, x_);
		for (size_t i = 0; i < m_number_variables; ++i) {
			x[i] = 0;
			for (size_t j = 0; j < m_number_variables; ++j) {
				x[i] += m_rotation[i][j] * x_[j];
			}
		}
		delete[] x_;
		x_ = 0;
	}
}