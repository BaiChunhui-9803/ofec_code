#include "dtlz.h"
#include <fstream>
#include "../../../../../core/global.h"

namespace ofec {
	void DTLZ::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeVariable(v.get<int>("number of variables"));//recommend Xm=5,n=M+Xm-1
		setDomain(0., 1.);

		resizeObjective(v.get<int>("number of objectives"));
		if (m_number_variables < m_number_objectives)
			throw MyExcept("The number of variables should be no less than the number of objectives.");
		for (size_t i = 0; i < m_number_objectives; ++i)
			m_optimize_mode[i] = OptimizeMode::kMinimize;

		loadParetoFront();
	}

	void DTLZ::generateParetoFront() {
		const std::string problem_name[] = { "MOP_DTLZ1", "MOP_DTLZ2", "MOP_DTLZ3", "MOP_DTLZ4", "MOP_DTLZ7" };
		const int M[5] = { 3, 5, 8, 10, 15 };
		bool flag1(false);
		for (int i = 0; i < 5; ++i) {
			if (m_name == problem_name[i]) {
				flag1 = true;
				break;
			}
		}
		if (!flag1)
			throw MyExcept("The problem name should be included in MOP_DTLZ1-MOP_DTLZ4, MOP_DTLZ7");
		else {
			bool flag2(false);
			for (int i = 0; i < 5; ++i) {
				if (m_number_objectives == M[i]) {
					flag2 = true;
					break;
				}
			}
			if (!flag2)
				throw MyExcept("The number of objectives should be included in 3,5,8,10,15");
			else { // generate PF
				std::stringstream os;
				os <<g_working_dir<< "/instance/problem/continuous/multi_objective/dtlz/data/PF_" << m_name << "_" << m_number_objectives << "objs.dtlz";
				std::ofstream ofile(os.str());

				if (m_number_objectives == 3)
					generateOnelayerPF(ofile, m_name, m_number_objectives, 12);
				else if (m_number_objectives == 5)
					generateOnelayerPF(ofile, m_name, m_number_objectives, 6);
				else if(m_number_objectives == 8 || m_number_objectives == 10)
					generateTwolayersPF(ofile, m_name, m_number_objectives, 3, 2);
				else
					generateTwolayersPF(ofile, m_name, m_number_objectives, 2, 1);
				ofile.close();
			}
		}

	}

	void DTLZ::loadParetoFront() {
		m_optima.reset(new Optima<>);
		std::stringstream os;
		os << g_working_dir<<"/instance/problem/continuous/multi_objective/dtlz/data/PF_" << m_name << "_" << m_number_objectives << "objs.dtlz";
		std::ifstream infile(os.str());
		if (!infile) {
			generateParetoFront();
			infile.close();
			infile.clear();
			infile.open(os.str());
		}
		std::string str;
		size_t line = 0;
		while (getline(infile, str))
			++line;
		m_optima->resizeObjectiveSet(line);
		infile.close();
		infile.clear();
		infile.open(os.str());
		for (size_t i = 0; i < line; i++) {
			std::vector<Real> temp_obj(m_number_objectives);
			for (size_t j = 0; j < m_number_objectives; j++)
				infile >> temp_obj[j];
			m_optima->setObjective(temp_obj, i);
		}
		m_optima->setObjectiveGiven(true);
		infile.close();
	}

	void DTLZ::generateRecursive(TFront * pf, TObjVec * pt, size_t number_objectives, size_t left, size_t total, size_t element) {
		if (element == number_objectives - 1) {
			(*pt)[element] = (Real)left;
			pf->push_back(*pt);
		}
		else {
			for (size_t i = 0; i <= left; i += 1) {
				(*pt)[element] = (Real)i;
				generateRecursive(pf, pt, number_objectives, left - i, total, element + 1);
			}
		}
	}

	void DTLZ::generateWeight(TFront * pf, size_t M, size_t p) {
		TObjVec pt(M);
		generateRecursive(pf, &pt, M, p, p, 0);
	}

	void DTLZ::generateOnelayerPF(std::ostream & os, const std::string & problem_name, int M, int p) {
		TFront PF;
		int num_objectives = M, num_divisions = p;
		generateWeight(&PF, num_objectives, num_divisions);
		if (problem_name == "MOP_DTLZ1") {
			for (size_t i = 0; i<PF.size(); i += 1) {
				for (size_t j = 0; j<PF[i].size(); j += 1)
					os << (0.5*PF[i][j]) / num_divisions << ' ';
				os << std::endl;
			}
		}
		else { // DTLZ2-4 
			for (size_t i = 0; i<PF.size(); i += 1) {
				Real sum = 0;
				for (size_t j = 0; j<PF[i].size(); j += 1)
					sum += PF[i][j] * PF[i][j];
				Real k = sqrt(1.0 / sum);
				for (size_t j = 0; j<PF[i].size(); j += 1)
					os << k*PF[i][j] << ' ';
				os << std::endl;
			}
		}
	}

	void DTLZ::generateTwolayersPF(std::ostream & os, const std::string & problem_name, int M, int outside_p, int inside_p) {
		generateOnelayerPF(os, problem_name, M, outside_p);
		TFront PF;
		int num_objectives = M, num_divisions = inside_p;
		generateWeight(&PF, num_objectives, num_divisions);
		for (size_t i = 0; i<PF.size(); i += 1) {
			for (size_t j = 0; j<PF[i].size(); j += 1)
				PF[i][j] = (static_cast<Real>(num_divisions) / M + PF[i][j]) / 2; // (k=num_divisions/M, k, k, ..., k) is the center point
		}
		if (problem_name == "DTLZ1") {
			for (size_t i = 0; i<PF.size(); i += 1) {
				for (size_t j = 0; j<PF[i].size(); j += 1)
					os << (0.5*PF[i][j]) / num_divisions << ' ';
				os << std::endl;
			}
		}
		else {// DTLZ2-4
			for (size_t i = 0; i < PF.size(); i += 1) {
				Real sum = 0;
				for (size_t j = 0; j < PF[i].size(); j += 1)
					sum += PF[i][j] * PF[i][j];
				Real k = sqrt(1.0 / sum);
				for (size_t j = 0; j < PF[i].size(); j += 1)
					os << k*PF[i][j] << ' ';
				os << std::endl;
			}
		}
	}
}
