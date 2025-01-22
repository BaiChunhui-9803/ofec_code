#include "zdt.h"
#include <fstream>
#include "../../../../../core/global.h"
#include"../../../../../utility/nondominated_sorting/fast_sort.h"

namespace ofec {
	void ZDT::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;;
		resizeObjective(v.get<int>("number of objectives"));
		m_num_reference_points = v.get<int>("# reference points");
		if (m_number_objectives != 2)
			throw MyExcept("The number of objectives must be equal to 2");
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		m_optimize_mode[1] = OptimizeMode::kMinimize;

		resizeVariable(v.get<int>("number of variables"));//recomend n=10
		if (m_number_variables < 2)
			throw MyExcept("The number of variables must be no less than 2");
		if (m_name == "MOP_ZDT4") {
			m_domain.setRange(0., 1., 0);
			for (size_t i = 1; i < m_number_variables; ++i) {
				m_domain.setRange(-5., 5., i);
			}
		}
		else {
			setDomain(0., 1.);
		}
		//m_domain_update = true;
	}

	//void ZDT1::sampleParetoFront(size_t sample_num) {
	//	size_t num = sample_num;
	//	std::vector<std::vector<Real>> all_objs;
	//	if (m_name == "MOP_ZDT3" || m_name == "MOP_ZDT6") {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			if (m_name == "MOP_ZDT3") {
	//				for (size_t j = 0; j < m_number_objectives; j++) {
	//					if (j == 0) {
	//						temp_obj[j] = (Real)i / (num - 1.);
	//					}
	//					else {
	//						temp_obj[j] = 1 - std::sqrt(temp_obj[0]) - temp_obj[0] * std::sin(10 * OFEC_PI * temp_obj[0]);
	//					}
	//				}
	//			}
	//			else if (m_name == "MOP_ZDT6") {
	//				Real temp = (Real)i / (num - 1.);
	//				for (size_t j = 0; j < m_number_objectives; j++) {
	//					if (j == 0) {
	//						temp_obj[j] = 1 - std::exp(-4 * temp) * std::pow(std::sin(6 * OFEC_PI * temp), 6);
	//					}
	//					else {
	//						temp_obj[j] = 1 - std::pow(temp_obj[0], 2);
	//					}
	//				}
	//			}
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
	//			if (rank[i] == 0) {
	//				m_optima->appendObj(all_objs[i]);
	//			}
	//		}
	//	}
	//	else {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			for (size_t j = 0; j < m_number_objectives; j++) {
	//				if (j == 0) {
	//					temp_obj[j] = (Real)i / (num - 1.);
	//				}
	//				else {
	//					if (m_name == "MOP_ZDT1" || m_name == "MOP_ZDT4") {
	//						temp_obj[j] = 1 - std::sqrt(temp_obj[0]);
	//					}
	//					else if (m_name == "MOP_ZDT2") {
	//						temp_obj[j] = 1 - std::pow(temp_obj[0], 2);
	//					}
	//				}
	//			}
	//			all_objs.emplace_back(temp_obj);
	//		}
	//		for (size_t i = 0; i < all_objs.size(); ++i) {
	//			m_optima->appendObj(all_objs[i]);
	//		}
	//	}
	//	m_optima->setObjectiveGiven(true);

	//	//write into file
	//	std::ofstream out1;
	//	std::stringstream os1;
	//	os1 << g_working_dir << "/instance/problem/continuous/multi_objective/zdt/data/" << m_name << "_obj.txt";
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

	//void ZDT1::loadParetoFront(size_t sample_num) {
	//	std::stringstream os1;
	//	os1 << g_working_dir << "/instance/problem/continuous/multi_objective/zdt/data/" << m_name << "_obj.txt";
	//	auto ss = os1.str();
	//	std::ifstream infile(ss);
	//	if (!infile.is_open())
	//	{
	//		throw MyExcept("open PF file of zdt problem is fail");
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
	//			throw MyExcept("open PF file of zdt problem is fail");
	//		}
	//		//m_optima.reset(new Optima<>());
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

	//void ZDT1::sampleParetoSets(size_t sample_num) {
	//	size_t num = sample_num;
	//	std::vector<VariableVector<Real>> all_vars;
	//	std::vector<std::vector<Real>> all_objs;
	//	if (m_name == "MOP_ZDT3" || m_name == "MOP_ZDT6") {
	//		for (size_t i = 0; i < num; i++) {
	//			VariableVector<Real> temp_var(m_number_variables);
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			for (size_t j = 0; j < m_number_variables; j++) {
	//				if (j == 0) {
	//					temp_var[j] = (Real)i / (num - 1.);
	//				}
	//				else {
	//					temp_var[j] = 0.;
	//				}
	//			}
	//			evaluateObjective(temp_var.data(), temp_obj);
	//			all_vars.emplace_back(temp_var);
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
	//			if (rank[i] == 0) {
	//				dynamic_cast<Optima<>&>(*m_optima).appendVar(all_vars[i]);
	//			}
	//		}
	//	}
	//	else {
	//		for (size_t i = 0; i < num; i++) {
	//			VariableVector<Real> temp_var(m_number_variables);
	//			for (size_t j = 0; j < m_number_variables; j++) {
	//				if (j == 0) {
	//					temp_var[j] = (Real)i / (num - 1.);
	//				}
	//				else {
	//					temp_var[j] = 0.;
	//				}
	//			}
	//			all_vars.emplace_back(temp_var);
	//		}
	//		for (size_t i = 0; i < all_vars.size(); ++i) {
	//			dynamic_cast<Optima<>&>(*m_optima).appendVar(all_vars[i]);
	//		}
	//	}
	//	m_optima->setVariableGiven(true);
	//}
}