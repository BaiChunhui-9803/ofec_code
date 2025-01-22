#include "moeadde_f.h"
#include "../../../../../core/global.h"
#include <fstream>

namespace ofec {
	void MOEADDE_F::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeVariable(v.get<int>("number of variables"));
		m_num_reference_points = v.get<int>("# reference points");
		m_domain.setRange(0., 1., 0);
		if (m_name == "MOP_MOEADDE_F1" || m_name == "MOP_MOEADDE_F7" || m_name == "MOP_MOEADDE_F8") {
			for (size_t i = 1; i < m_number_variables; i++)
				m_domain.setRange(0., 1., i);
		}
		if (m_name == "MOP_MOEADDE_F2" || m_name == "MOP_MOEADDE_F3"|| m_name == "MOP_MOEADDE_F4" || m_name == "MOP_MOEADDE_F5"|| m_name == "MOP_MOEADDE_F9") {
			for (size_t i = 1; i < m_number_variables; i++)
				m_domain.setRange(-1., 1., i);
		}
		else if(m_name == "MOP_MOEADDE_F6"){
			m_domain.setRange(0., 1., 1);
			for (size_t i = 2; i < m_number_variables; i++)
				m_domain.setRange(-2., 2., i);
		}

		//m_domain_update = true;

		resizeObjective(v.get<int>("number of objectives"));
		if (m_name == "MOP_MOEADDE_F6") {
			if (m_number_objectives != 3)
				throw MyExcept("The number of objectives must be equal to 3.");
		}
		else {
			if (m_number_objectives != 2)
				throw MyExcept("The number of objectives must be equal to 2.");
		}
		for (size_t i = 0; i < m_number_objectives; ++i)
			m_optimize_mode[i] = OptimizeMode::kMinimize;

		//loadParetoSets();
	}

	//void MOEADDE_F::loadParetoFront(size_t sample_num) {
	//	std::stringstream os1;
	//	os1 << g_working_dir << "/instance/problem/continuous/multi_objective/moeadde_f/data/" << m_name << "_obj.txt";
	//	auto ss = os1.str();
	//	std::ifstream infile(ss);
	//	if (!infile.is_open())
	//	{
	//		throw MyExcept("open PF file of MOEADDE_F problem is fail");
	//	}
	//	std::string str;
	//	size_t line = 0;
	//	while (getline(infile, str))
	//		++line;
	//	infile.close();
	//	infile.clear();
	//	if (line == 0) {
	//		sampleParetoFront(10000);
	//	}
	//	else {
	//		infile.open(os1.str());
	//		if (!infile) {
	//			throw MyExcept("open PF file of MOEADDE_F problem is fail");
	//		}
	//		m_optima->resizeObjectiveSet(line);
	//		for (size_t i = 0; i < line; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			for (size_t j = 0; j < m_number_objectives; j++)
	//				infile >> temp_obj[j];
	//			m_optima->setObjective(temp_obj, i);
	//		}
	//		m_optima->setObjectiveGiven(true);
	//		infile.close();
	//		infile.clear();
	//	}
	//}

	//void MOEADDE_F::sampleParetoFront(size_t sample_num) {
	//	size_t num = sample_num;
	//	std::vector<std::vector<Real>> all_objs;
	//	if (m_name == "MOP_MOEADDE_F9") {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			for (size_t j = 0; j < m_number_objectives; j++) {
	//				if (j == 0)
	//					temp_obj[j] = (Real)i / (num - 1.);
	//				else
	//					temp_obj[j] = 1 - temp_obj[0]*temp_obj[0];
	//			}
	//			all_objs.emplace_back(temp_obj);
	//		}
	//	}
	//	else if (m_name == "MOP_MOEADDE_F6") {
	//		size_t div_num = 100;
	//		for (size_t i = 0; i < div_num; i++) {
	//			Real temp1 = (Real)i / (div_num - 1.);
	//			for (size_t j = 0; j < div_num; j++) {
	//				Real temp2 = (Real)j / (div_num - 1.);
	//				std::vector<Real> temp_obj(m_number_objectives);
	//				temp_obj[0] = std::cos(temp1 * OFEC_PI / 2) * std::cos(temp2 * OFEC_PI / 2);
	//				temp_obj[1] = std::cos(temp1 * OFEC_PI / 2) * std::sin(temp2 * OFEC_PI / 2);
	//				temp_obj[2] = std::sin(temp1 * OFEC_PI / 2);
	//				all_objs.emplace_back(temp_obj);
	//			}
	//		}
	//	}
	//	else {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			temp_obj[0] = (Real)i / (num - 1.);
	//			temp_obj[1] = 1 - std::sqrt(temp_obj[0]);
	//			all_objs.emplace_back(temp_obj);
	//		}
	//	}
	//	for (size_t i = 0; i < all_objs.size(); ++i) {
	//		m_optima->appendObj(all_objs[i]);
	//	}
	//	m_optima->setObjectiveGiven(true);

	//	//write into file
	//	std::ofstream out1;
	//	std::stringstream os1;
	//	os1 << g_working_dir << "/instance/problem/continuous/multi_objective/moeadde_f/data/" << m_name << "_obj.txt";
	//	out1.open(os1.str());
	//	if (!out1) {
	//		throw("can not open file");
	//	}
	//	for (size_t i = 0; i < m_optima->numberObjectives(); ++i) {
	//		auto& objs = m_optima->objective(i);
	//		for (size_t j = 0; j < objs.size(); ++j) {
	//			out1 << objs[j] << " ";
	//		}
	//		out1 << std::endl;
	//	}
	//	out1.close();
	//}

	std::vector<Real> MOEADDE_F::createVar(const std::vector<Real>& s) const{
		//根据不同的问题，构造一个最优解
		//auto curSol = new Solution<>(numberObjectives(), numberConstraints(), numberVariables());
		//new_var = dynamic_cast<const VariableVector<Real>&>(s);
		std::vector<Real> new_var;
		for (size_t i = 0; i < s.size(); ++i) {
			new_var.push_back(s[i]);
		}
		
		if (m_name != "MOP_MOEADDE_F6") {
			for (size_t j = 1; j < m_number_variables; j++) {
				if (m_name == "MOP_MOEADDE_F1" || m_name == "MOP_MOEADDE_F7" || m_name == "MOP_MOEADDE_F8") {
					if (m_number_variables == 2)
						new_var[j] = std::pow(new_var[0], 0.5 * (1. + 3 * j / (m_number_variables - 1)));
					else
						new_var[j] = std::pow(new_var[0], 0.5 * (1. + 3 * (j - 1) / (m_number_variables - 2)));
				}
				else if (m_name == "MOP_MOEADDE_F2") {
					new_var[j] = std::sin(4 * OFEC_PI * new_var[0] + (j + 1) * OFEC_PI / m_number_variables);
				}
				else if (m_name == "MOP_MOEADDE_F3") {
					if (j % 2 == 1)
						new_var[j] = 0.8 * new_var[0] * std::sin(6 * OFEC_PI * new_var[0] + (j + 1) * OFEC_PI / m_number_variables);
					else
						new_var[j] = 0.8 * new_var[0] * std::cos(6 * OFEC_PI * new_var[0] + (j + 1) * OFEC_PI / m_number_variables);
				}
				else if (m_name == "MOP_MOEADDE_F4") {
					if (j % 2 == 1)
						new_var[j] = 0.8 * new_var[0] * std::sin(6 * OFEC_PI * new_var[0] + (j + 1) * OFEC_PI / m_number_variables);
					else
						new_var[j] = 0.8 * new_var[0] * std::cos((6 * OFEC_PI * new_var[0] + (j + 1) * OFEC_PI / m_number_variables) / 3);
				}
				else if (m_name == "MOP_MOEADDE_F5") {
					if (j % 2 == 1)
						new_var[j] = (0.3 * std::pow(new_var[0], 2) * std::cos(24 * OFEC_PI * new_var[0] + 4 * (j + 1) * OFEC_PI / m_number_variables) + 0.6 * new_var[0]) * std::sin(6 * OFEC_PI * new_var[0] + (j + 1) * OFEC_PI / m_number_variables);
					else
						new_var[j] = (0.3 * std::pow(new_var[0], 2) * std::cos(24 * OFEC_PI * new_var[0] + 4 * (j + 1) * OFEC_PI / m_number_variables) + 0.6 * new_var[0]) * std::cos(6 * OFEC_PI * new_var[0] + (j + 1) * OFEC_PI / m_number_variables);
				}
				else if (m_name == "MOP_MOEADDE_F9") {
					new_var[j] = std::sin(6 * OFEC_PI * new_var[0] + (j + 1) * OFEC_PI / m_number_variables);
				}
			}
		}
		else {
			for (size_t k = 2; k < m_number_variables; k++) {
				new_var[k] = 2 * new_var[1] * std::sin(2 * OFEC_PI * new_var[0] + (1 + k) * OFEC_PI / m_number_variables);
			}
		}
		return new_var;
	}
}