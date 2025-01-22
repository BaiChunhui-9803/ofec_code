#include "moeadm2m_f.h"
#include "../../../../../core/global.h"
#include <fstream>

namespace ofec {
	void MOEADM2M_F::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeVariable(v.get<int>("number of variables"));
		m_num_reference_points = v.get<int>("# reference points");
		for (size_t i = 0; i < m_number_variables; i++) {
			m_domain.setRange(0., 1., i);
		}

		//m_domain_update = true;

		resizeObjective(v.get<int>("number of objectives"));
		if (m_name == "MOP_MOEADM2M_F6" || m_name == "MOP_MOEADM2M_F7") {
			if (m_number_objectives != 3)
				throw MyExcept("The number of objectives must be equal to 3.");
			if (m_number_variables < 2)
				throw MyExcept("The number of variables must be over 1.");
		}
		else {
			if (m_number_objectives != 2)
				throw MyExcept("The number of objectives must be equal to 2.");
		}
		for (size_t i = 0; i < m_number_objectives; ++i)
			m_optimize_mode[i] = OptimizeMode::kMinimize;

		//loadParetoSets();
	}

	//void MOEADM2M_F::loadParetoFront(size_t sample_num) {
	//	std::stringstream os1;
	//	os1 << g_working_dir << "/instance/problem/continuous/multi_objective/moeadm2m_f/data/" << m_name << "_obj.txt";
	//	auto ss = os1.str();
	//	std::ifstream infile(ss);
	//	if (!infile.is_open())
	//	{
	//		throw MyExcept("open PF file of MOEADM2M_F problem is fail");
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
	//			throw MyExcept("open PF file of MOEADM2M_F problem is fail");
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

	//void MOEADM2M_F::sampleParetoFront(size_t sample_num) {
	//	size_t num = sample_num;
	//	std::vector<std::vector<Real>> all_objs;
	//	if (m_name == "MOP_MOEADM2M_F1" || m_name == "MOP_MOEADM2M_F5") {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			temp_obj[0] = (Real)i / (num - 1.);
	//			temp_obj[1] = 1 - std::sqrt(temp_obj[0]);
	//			all_objs.emplace_back(temp_obj);
	//		}
	//		for (size_t i = 0; i < all_objs.size(); ++i) {
	//			m_optima->appendObj(all_objs[i]);
	//		}
	//	}
	//	else if (m_name == "MOP_MOEADM2M_F2"){
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			temp_obj[0] = (Real)i / (num - 1.);
	//			temp_obj[1] = 1 - std::pow(temp_obj[0],2);
	//			all_objs.emplace_back(temp_obj);
	//		}
	//		for (size_t i = 0; i < all_objs.size(); ++i) {
	//			m_optima->appendObj(all_objs[i]);
	//		}
	//	}
	//	else if (m_name == "MOP_MOEADM2M_F3") {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			temp_obj[0] = (Real)i / (num - 1.);
	//			temp_obj[1] = std::sqrt(1 - std::pow(temp_obj[0], 2));
	//			all_objs.emplace_back(temp_obj);
	//		}
	//		for (size_t i = 0; i < all_objs.size(); ++i) {
	//			m_optima->appendObj(all_objs[i]);
	//		}
	//	}
	//	else if (m_name == "MOP_MOEADM2M_F4") {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			temp_obj[0] = (Real)i / (num - 1.);
	//			temp_obj[1] = 1-std::pow(temp_obj[0],0.5) * std::pow(std::cos(2*OFEC_PI*temp_obj[0]), 2);
	//			all_objs.emplace_back(temp_obj);
	//		}
	//		//first nd-sort, then add into optima
	//		std::vector<std::vector<Real>*> objs;
	//		for (size_t i = 0; i < all_objs.size(); ++i) {
	//			objs.emplace_back(&all_objs[i]);
	//		}
	//		std::vector<int> rank;
	//		ofec::nd_sort::fastSort<Real>(objs, rank, m_optimize_mode);
	//		for (size_t i = 0; i < rank.size(); ++i) {
	//			if (rank[i] == 0)
	//				m_optima->appendObj(all_objs[i]);
	//		}
	//	}
	//	else if (m_name == "MOP_MOEADM2M_F6") {
	//		size_t div_num = 100;
	//		for (size_t i = 0; i < div_num; i++) {
	//			Real temp1 = (Real)i / (div_num - 1.);
	//			for (size_t j = 0; j < div_num; j++) {
	//				Real temp2 = (Real)j / (div_num - 1.);
	//				std::vector<Real> temp_obj(m_number_objectives);
	//				temp_obj[0] = temp1 * temp2;
	//				temp_obj[1] = temp1 * (1-temp2);
	//				temp_obj[2] = 1 - temp1;
	//				all_objs.emplace_back(temp_obj);
	//			}
	//		}
	//		for (size_t i = 0; i < all_objs.size(); ++i) {
	//			m_optima->appendObj(all_objs[i]);
	//		}
	//	}
	//	else if (m_name == "MOP_MOEADM2M_F7") {
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
	//		for (size_t i = 0; i < all_objs.size(); ++i) {
	//			m_optima->appendObj(all_objs[i]);
	//		}
	//	}

	//	m_optima->setObjectiveGiven(true);

	//	//write into file
	//	std::ofstream out1;
	//	std::stringstream os1;
	//	os1 << g_working_dir << "/instance/problem/continuous/multi_objective/moeadm2m_f/data/" << m_name << "_obj.txt";
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

	std::vector<Real> MOEADM2M_F::createVar(const std::vector<Real>& s) const {
		//根据不同的问题，构造一个最优解
		//auto curSol = new Solution<>(numberObjectives(), numberConstraints(), numberVariables());
		//new_var = dynamic_cast<const VariableVector<Real>&>(s);
		std::vector<Real> new_var;
		for (size_t i = 0; i < s.size(); ++i) {
			new_var.push_back(s[i]);
		}

		if (m_name == "MOP_MOEADM2M_F6" || m_name == "MOP_MOEADM2M_F7") {
			for (size_t k = 2; k < m_number_variables; k++) {
				new_var[k] = new_var[0] * new_var[1];
			}
		}
		else {
			for (size_t j = 1; j < m_number_variables; j++) {
				new_var[j] = std::sin(0.5 * OFEC_PI * new_var[0]);
			}
		}
		return new_var;
	}

//	void MOEADM2M_F::updateOptima() {
//		m_optima.reset(new Optima<>);
//		loadParetoFront(10000);
//#ifdef OFEC_DEMO
//		sampleParetoSets(10000);
//#endif //
//	}
//
//	void MOEADM2M_F::sampleParetoSols(size_t sample_num) {
//		size_t num = sample_num;
//		std::vector<VariableVector<Real>> all_vars;
//		std::vector<std::vector<Real>> all_objs;
//		if (m_name == "MOP_MOEADM2M_F6" || m_name == "MOP_MOEADM2M_F7") {
//			size_t div_num = 100;
//			for (size_t i = 0; i < div_num; i++) {
//				for (size_t j = 0; j < div_num; j++) {
//					VariableVector<Real> temp_var(m_number_variables);
//					for (size_t k = 0; k < m_number_variables; k++) {
//						if (k == 0) {
//							temp_var[k] = (Real)i / (div_num - 1.);
//						}
//						else if (k == 1) {
//							temp_var[k] = (Real)j / (div_num - 1.);
//						}
//						else {
//							temp_var[k] = temp_var[0] * temp_var[1];
//						}
//					}
//					all_vars.emplace_back(temp_var);
//				}
//			}
//		}
//		else {
//			for (size_t i = 0; i < num; i++) {
//				VariableVector<Real> temp_var(m_number_variables);
//				for (size_t j = 0; j < m_number_variables; j++) {
//					if (j == 0)
//						temp_var[j] = (Real)i / (num - 1.);
//					else {
//						temp_var[j] = std::sin(0.5 * OFEC_PI * temp_var[0]);
//					}
//				}
//				all_vars.emplace_back(temp_var);
//			}
//		}
//		for (size_t i = 0; i < all_vars.size(); ++i) {
//			dynamic_cast<Optima<>&>(*m_optima).appendVar(all_vars[i]);
//		}
//		m_optima->setVariableGiven(true);
//	}
}