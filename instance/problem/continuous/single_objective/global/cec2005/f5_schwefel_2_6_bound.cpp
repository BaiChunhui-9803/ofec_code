#include "f5_schwefel_2_6_bound.h"
#include "../../../../../../core/global.h"
#include <sstream>

namespace ofec::cec2005 {
	void Schwefel_2_6_Bound::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void Schwefel_2_6_Bound::initialize_(Environment *env) {
		Function::initialize_(env);
		resizeVariable(m_number_variables);
		setDomain(-100, 100);
		loadData("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005_F05"); // load data into m_a and m_b
		setBias(-310);
		setGlobalOpt(m_translation.data());
		m_objective_accuracy = 1.0e-6;
	}

	void Schwefel_2_6_Bound::loadData(const std::string &path) {
		std::string sa;
		std::stringstream ss;
		ss << m_number_variables << "Dim.txt";
		sa = ss.str();
		sa.insert(0, "_a_");
		sa.insert(0, path);
		sa.insert(0, g_working_directory);// data path

		m_a.assign(m_number_variables, std::vector<Real>(m_number_variables));
		m_b.resize(m_number_variables);

		std::ifstream in_a(sa);
		if (in_a.fail()) {
			for (int i = 0; i < m_number_variables; ++i) {
				for (int j = 0; j < m_number_variables; ++j) {
					m_a[i][j] = int(-500.0 + m_random->uniform.next() * 1000);
				}
			}
			std::ofstream out(sa.c_str());
			for (int i = 0; i < m_number_variables; ++i) {
				for (int j = 0; j < m_number_variables; j++) {
					out << m_a[i][j] << " ";
				}
			}
			out.close();
		}
		else {
			std::string row;
			for (int i = 0; i < m_number_variables; ++i) {
				std::getline(in_a, row);
				std::stringstream sstr_row(row);
				for (int j = 0; j < m_number_variables; j++) {
					sstr_row >> m_a[i][j];
				}
			}
		}
		in_a.close();

		setOriginalGlobalOpt();


		m_translation.resize(m_number_variables);
		std::string so;
		so = ss.str();
		so.insert(0, "_Opt_");

		so.insert(0, path);
		so.insert(0, g_working_directory);// data path

		std::ifstream in;
		in.open(so.data());
		if (in.fail()) {
			setTranslation(m_original_optima.solution(0).variable().data(), m_random.get());
			for (int i = 0; i < m_number_variables; ++i) {
				if (i < m_number_variables / 4) m_translation[i] = -100;
				else if (i >= m_number_variables * 3 / 4 - 1) m_translation[i] = 100;
			}

			std::ofstream out(so.c_str());
			for (int j = 0; j < m_number_variables; j++)        out << m_translation[j] << " ";
			out.close();
		}
		else {
			for (int j = 0; j < m_number_variables; j++) {
				in >> m_translation[j];
			}
		}
		in.close();

		for (size_t i = 0; i < m_number_variables / 4 + 1; ++i)
			m_translation[i] = -100;
		for (size_t i = m_number_variables * 3 / 4; i < m_number_variables; ++i)
			m_translation[i] = 100;

		for (int i = 0; i < m_number_variables; ++i) {
			m_b[i] = 0;
			for (int j = 0; j < m_number_variables; ++j) {
				m_b[i] += m_a[i][j] * m_translation[j];
			}
		}
	}

	void Schwefel_2_6_Bound::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		std::vector<Real> temp_vector(m_number_variables);
		for (int i = 0; i < m_number_variables; ++i) {
			for (int j = 0; j < m_number_variables; ++j) {
				temp_vector[j] = 0;
				for (int k = 0; k < m_number_variables; ++k) {
					temp_vector[j] += m_a[j][k] * x[k];
				}
			}
			for (int j = 0; j < m_number_variables; ++j) {
				temp_vector[j] -= m_b[j];
			}
		}
		Real temp_max = abs(temp_vector[0]);
		for (size_t j = 1; j < m_number_variables; ++j) {
			if (abs(temp_vector[j]) > temp_max)
				temp_max = abs(temp_vector[j]);
		}
		obj[0] = temp_max + m_bias;
	}
}