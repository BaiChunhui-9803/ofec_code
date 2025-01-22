#include "zdt1.h"

namespace ofec {
	void ZDT1::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Real g = 0;
		for (size_t n = 1; n < m_number_variables; n++) {
			g = g + x[n];
		}
		g = 1 + 9 * g / (m_number_variables - 1);
		obj[0] = x[0];
		obj[1] = g * (1 - pow(obj[0] / g, 0.5));
	}

	void ZDT1::updateOptima() {
		m_optima.reset(new Optima<>);
		//loadParetoFront(10000);
		sampleParetoSols(m_num_reference_points);
//#ifdef OFEC_DEMO
//		sampleParetoSets(10000);
//#endif //
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

	void ZDT1::sampleParetoSols(size_t sample_num) {
		std::vector<VariableVector<Real>> all_vars;
		std::vector<std::vector<Real>> all_objs;
		for (size_t i = 0; i < sample_num; i++) {
			VariableVector<Real> temp_var(m_number_variables);
			std::vector<Real> temp_obj(m_number_objectives);
			for (size_t j = 0; j < m_number_variables; j++) {
				if (j == 0) {
					temp_var[j] = (Real)i / (sample_num - 1.);
				}
				else {
					temp_var[j] = 0.;
				}
			}
			Solution<> sol(m_number_objectives, m_number_constraints);
			sol.variable() = temp_var;
			evaluateObjective(temp_var.data(), temp_obj);
			sol.objective() = temp_obj;
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(sol);
		}
	}
}