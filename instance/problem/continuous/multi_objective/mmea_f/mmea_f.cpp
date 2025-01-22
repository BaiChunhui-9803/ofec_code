#include "mmea_f.h"
#include "../../../../../core/global.h"
#include "../../../../../core/problem/solution.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include <fstream>

namespace ofec {
	void MMEA_F::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;

		resizeVariable(v.get<int>("number of variables"));
		m_num_reference_points = v.get<int>("# reference points");
		for (size_t i = 0; i < m_number_variables; i++) {
			m_domain.setRange(0., 1., i);
		}
			
		//m_domain_update = true;

		resizeObjective(v.get<int>("number of objectives"));
		if (m_name == "MOP_MMEA_F7") {
			if (m_number_objectives != 3)
				throw MyExcept("The number of objectives must be equal to 3.");
			if (m_number_variables < 3)
				throw MyExcept("The number of variables must be over 2.");
		}
		else {
			if (m_number_objectives != 2)
				throw MyExcept("The number of objectives must be equal to 2.");
			if (m_number_variables < 2)
				throw MyExcept("The number of variables must be over 1.");
			if (m_name == "MOP_MMEA_F6") {
				if (m_number_variables < 3)
					throw MyExcept("The number of variables must be over 2.");
			}
		}

		for (size_t i = 0; i < m_number_objectives; ++i)
			m_optimize_mode[i] = OptimizeMode::kMinimize;

		//loadParetoSets();
	}

	//void MMEA_F::loadParetoFront(size_t sample_num) {
	//	std::stringstream os1;
	//	os1 << g_working_dir << "/instance/problem/continuous/multi_objective/mmea_f/data/" << m_name << "_obj.txt";
	//	auto ss = os1.str();
	//	std::ifstream infile(ss);
	//	if (!infile.is_open())
	//	{
	//		throw MyExcept("open PF file of MMEA_F problem is fail");
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
	//			throw MyExcept("open PF file of MMEA_F problem is fail");
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

	//void MMEA_F::sampleParetoFront(size_t sample_num) {
	//	size_t num = sample_num;
	//	std::vector<std::vector<Real>> all_objs;
	//	if (m_name == "MOP_MMEA_F3") {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			temp_obj[0] = (Real)i / (num - 1.);
	//			temp_obj[1] = 1 - std::sqrt(temp_obj[0]);
	//			all_objs.emplace_back(temp_obj);
	//		}
	//	}
	//	else if (m_name == "MOP_MMEA_F4" || m_name == "MOP_MMEA_F6") {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			temp_obj[0] = (Real)i / (num - 1.);
	//			temp_obj[1] = 1 - std::pow(temp_obj[0], 2);
	//			all_objs.emplace_back(temp_obj);
	//		}
	//	}
	//	else if (m_name == "MOP_MMEA_F5") {
	//		for (size_t i = 0; i < num; i++) {
	//			std::vector<Real> temp_obj(m_number_objectives);
	//			temp_obj[0] = (Real)i / (num - 1.);
	//			temp_obj[1] = 1 - temp_obj[0] + std::sin(2 * OFEC_PI * temp_obj[0])/2./OFEC_PI;
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
	//	else {
	//		size_t div_num = 25;
	//		for (size_t i = 0; i < div_num; i++) {
	//			Real temp1 = (Real)i / (div_num - 1.);
	//			for (size_t j = 0; j < div_num; j++) {
	//				Real temp2 = (Real)j / (div_num - 1.);
	//				for (size_t k = 0; k < div_num; k++) {
	//					Real temp3 = (Real)k / (div_num - 1.);
	//					std::vector<Real> temp_obj(m_number_objectives);
	//					temp_obj[0] = std::cos(0.25 * OFEC_PI *(temp1+temp2)) * std::sin(0.5 * OFEC_PI *temp3);
	//					temp_obj[1] = std::cos(0.25 * OFEC_PI * (temp1 + temp2)) * std::cos(0.5 * OFEC_PI * temp3);
	//					temp_obj[2] = std::sin(0.25 * OFEC_PI * (temp1 + temp2));
	//					all_objs.emplace_back(temp_obj);
	//				}
	//			}
	//		}
	//	}
	//	for (size_t i = 0; i < all_objs.size(); ++i) {
	//		m_optima->appendObj(all_objs[i]);
	//	}
	//	m_optima->setObjectiveGiven(true);

	//	//write into file
	//	std::ofstream out1;
	//	std::stringstream os1;
	//	os1 << g_working_dir << "/instance/problem/continuous/multi_objective/mmea_f/data/" << m_name << "_obj.txt";
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

	std::vector<Real> MMEA_F::createVar(const std::vector<Real>& s) const {
		//根据不同的问题，构造一个最优解
		//auto curSol = new Solution<>(numberObjectives(), numberConstraints(), numberVariables());
		//new_var = dynamic_cast<const VariableVector<Real>&>(s);
		std::vector<Real> temp_var;
		for (size_t i = 0; i < s.size(); ++i) {
			temp_var.push_back(s[i]);
		}
		if (m_name != "MOP_MMEA_F7" && m_name != "MOP_MMEA_F6") {
			for (size_t j = 2; j < m_number_variables; j++) {
				Real y = (temp_var[0] + temp_var[1]) / 2.;
				if (m_name == "MOP_MMEA_F3") {
					if (j % 2) {
						temp_var[j] = 0.5 + 0.5 * std::sin(0.5 * OFEC_PI * y) * std::cos(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
					}
					else {
						temp_var[j] = 0.5 + 0.5 * std::cos(0.5 * OFEC_PI * y) * std::sin(1. / 3 * (2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables));
					}
				}
				else if (m_name == "MOP_MMEA_F4") {
					if (j % 2) {
						temp_var[j] = 0.5 + 0.5 * y * std::cos(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
					}
					else {
						temp_var[j] = 0.5 + 0.5 * y * std::sin(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
					}
				}
				else if (m_name == "MOP_MMEA_F5") {
					if (j % 2) {
						temp_var[j] = 0.5 + 0.5 * y * std::cos(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
					}
					else {
						temp_var[j] = 0.5 + 0.5 * y * std::sin(1. / 3 * (2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables));
					}
				}
			}
		}
		else {
			for (size_t p = 3; p < m_number_variables; p++) {
				Real y = (temp_var[0] + temp_var[1] + temp_var[2]) / 3.;
				if (p % 2) {
					temp_var[p] = 0.5 + 0.5 * std::sin(0.5 * OFEC_PI * y) * std::cos(2 * OFEC_PI * y + (p + 1) * OFEC_PI / m_number_variables);
				}
				else {
					temp_var[p] = 0.5 + 0.5 * std::cos(0.5 * OFEC_PI * y) * std::sin(1. / 3 * (2 * OFEC_PI * y + (p + 1) * OFEC_PI / m_number_variables));
				}
			}
		}

		return temp_var;
	}

//	void MMEA_F::updateOptima() {
//		m_optima.reset(new Optima<>);
//		loadParetoFront(10000);
//#ifdef OFEC_DEMO
//		sampleParetoSets(100);
//#endif //
//	}
//
//	void MMEA_F::sampleParetoSets(size_t sample_num) {
//		size_t div_num = sample_num;
//		std::vector<VariableVector<Real>> all_vars;
//		std::vector<std::vector<Real>> all_objs;
//		if (m_name != "MOP_MMEA_F7" && m_name != "MOP_MMEA_F6") {
//			for (size_t i = 0; i < div_num; i++) {
//				for (size_t k = 0; k < div_num; ++k) {
//					VariableVector<Real> temp_var(m_number_variables);
//					for (size_t j = 0; j < m_number_variables; j++) {
//						if (j == 0)
//							temp_var[j] = (Real)i / (div_num - 1.);
//						else if (j == 1) {
//							temp_var[j] = (Real)k / (div_num - 1.);
//						}
//						else {
//							Real y = (temp_var[0] + temp_var[1]) / 2.;
//							if (m_name == "MOP_MMEA_F3") {
//								if (j % 2) {
//									temp_var[j] = 0.5 + 0.5 * std::sin(0.5 * OFEC_PI * y) * std::cos(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
//								}
//								else {
//									temp_var[j] = 0.5 + 0.5 * std::cos(0.5 * OFEC_PI * y) * std::sin(1. / 3 * (2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables));
//								}
//							}
//							else if (m_name == "MOP_MMEA_F4") {
//								if (j % 2) {
//									temp_var[j] = 0.5 + 0.5 * y * std::cos(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
//								}
//								else {
//									temp_var[j] = 0.5 + 0.5 * y * std::sin(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
//								}
//							}
//							else if (m_name == "MOP_MMEA_F5") {
//								if (j % 2) {
//									temp_var[j] = 0.5 + 0.5 * y * std::cos(2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables);
//								}
//								else {
//									temp_var[j] = 0.5 + 0.5 * y * std::sin(1. / 3 * (2 * OFEC_PI * y + (j + 1) * OFEC_PI / m_number_variables));
//								}
//							}
//						}
//					}
//					all_vars.emplace_back(temp_var);
//				}
//			}
//		}
//		else {
//			size_t div_num = 25;
//			for (size_t i = 0; i < div_num; i++) {
//				for (size_t j = 0; j < div_num; j++) {
//					for (size_t k = 0; k < div_num; k++) {
//						VariableVector<Real> temp_var(m_number_variables);
//						for (size_t p = 0; p < m_number_variables; p++) {
//							if (p == 0) {
//								temp_var[p] = (Real)i / (div_num - 1.);
//							}
//							else if (p == 1) {
//								temp_var[p] = (Real)j / (div_num - 1.);
//							}
//							else if (p == 2) {
//								temp_var[p] = (Real)k / (div_num - 1.);
//							}
//							else {
//								Real y = (temp_var[0]+temp_var[1]+temp_var[2]) / 3.;
//								if (p % 2) {
//									temp_var[p] = 0.5 + 0.5 * std::sin(0.5 * OFEC_PI * y) * std::cos(2 * OFEC_PI * y + (p + 1) * OFEC_PI / m_number_variables);
//								}
//								else {
//									temp_var[p] = 0.5 + 0.5 * std::cos(0.5 * OFEC_PI * y) * std::sin(1. / 3 * (2 * OFEC_PI * y + (p + 1) * OFEC_PI / m_number_variables));
//								}
//							}
//							all_vars.emplace_back(temp_var);
//						}
//					}
//				}
//			}
//		}
//		for (size_t i = 0; i < all_vars.size(); ++i) {
//			dynamic_cast<Optima<>&>(*m_optima).appendVar(all_vars[i]);
//		}
//		m_optima->setVariableGiven(true);
//	}
}