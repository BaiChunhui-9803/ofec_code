#include "quadratic_assignment.h"
#include "../../../../core/problem/solution.h"
#include "../../../../core/global.h"
#include <fstream>
#include <cstring>
#include <cstdlib>

namespace ofec {
	void QuadraticAssignment::initialize_() {
		Problem::initialize_();
		auto &v = *m_param;
		m_file_name = v.get<std::string>("dataFile1");
		readProblem();
		m_domain.resize(m_number_variables);
		for (size_t i = 0; i < m_number_variables; ++i)
			m_domain.setRange(0, m_number_variables - 1, i);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
	}

	void QuadraticAssignment::updateOptima() {
		m_optima.reset(new Optima<VariableVector<int>>());
		readOptima();
	}

	void QuadraticAssignment::evaluate_(SolutionBase &s) {
		VariableVector<int> &x = dynamic_cast<Solution<VariableVector<int>> &>(s).variable();
		std::vector<Real> &obj = dynamic_cast<Solution<VariableVector<int>> &>(s).objective();
		for (size_t i = 0; i < m_number_objectives; i++)
			obj[i] = 0;
		int row, col;
		for (size_t n = 0; n < m_number_objectives; n++) {
			for (size_t i = 0; i < m_number_variables; i++) {
				row = x[i];
				for (size_t j = 0; j < m_number_variables; j++) {
					col = x[j];
					obj[n] += mvv_distance[i][j] * mvv_flow[row][col];
				}
			}
		}
	}

	bool QuadraticAssignment::isValid(const SolutionBase &s) const {
		if (!m_if_valid_check)
			return true;
		const VariableVector<int> &x = dynamic_cast<const Solution<VariableVector<int>> &>(s).variable();
		for (int i = 0; i < m_number_variables; i++)  //judge the range		
			if ((x[i]) < m_domain.range(i).limit.first || (x[i]) > m_domain.range(i).limit.second)
				return false;
		std::vector<int> flag(m_number_variables, 0);  //judge whether has the same gene
		int temp;
		for (int i = 0; i < m_number_variables; i++) {
			temp = x[i];
			flag[temp] = 1;
		}
		for (int i = 0; i < m_number_variables; i++)
			if (flag[i] == 0)
				return false;
		return true;
	}

	void QuadraticAssignment::readProblem() {
#if defined(linux) || defined(__linux) || defined(__linux__) || defined(__APPLE__)
#define	strtok_s strtok_r
#endif
		size_t i;
		std::string Line;
		char *Keyword = nullptr;
		const char *Delimiters = " ():=\n\t\r\f\v\xef\xbb\xbf";
		std::ostringstream oss;
		std::ifstream infile;
		oss << static_cast<std::string>(g_working_dir);
		oss << "/instance/problem/combination/quadratic_assignment/" << m_file_name << ".qap";
		infile.open(oss.str().c_str());
		if (!infile) {
			throw MyExcept("read Quadratic Assignment data error");
		}
		char *savePtr;
		while (getline(infile, Line)) {
			if (!(Keyword = strtok_s((char *)Line.c_str(), Delimiters, &savePtr)))
				continue;
			for (i = 0; i < strlen(Keyword); i++)
				Keyword[i] = toupper(Keyword[i]);
			if (!strcmp(Keyword, "DIM")) {
				char *token = strtok_s(nullptr, Delimiters, &savePtr);
				m_number_variables = atoi(token);
			}
			else if (!strcmp(Keyword, "FLOW")) {
				mvv_flow.resize(m_number_variables);
				for (size_t i = 0; i < m_number_variables; i++)
					mvv_flow[i].resize(m_number_variables);
				for (size_t n = 0; n < m_number_variables; n++)
					for (i = 0; i < m_number_variables; i++)
						infile >> mvv_flow[n][i];
			}
			else if (!strcmp(Keyword, "DISTANCE")) {
				mvv_distance.resize(m_number_variables);
				for (size_t i = 0; i < m_number_variables; i++)
					mvv_distance[i].resize(m_number_variables);
				for (size_t n = 0; n < m_number_variables; n++)
					for (i = 0; i < m_number_variables; i++)
						infile >> mvv_distance[n][i];
			}
		}
		infile.close();
		infile.clear();
	}

	void QuadraticAssignment::readOptima() {
#if defined(linux) || defined(__linux) || defined(__linux__) || defined(__APPLE__)
#define	strtok_s strtok_r
#endif
		size_t i;
		std::string Line;
		char *Keyword = nullptr;
		const char *Delimiters = " ():=\n\t\r\f\v\xef\xbb\xbf";
		std::ostringstream oss;
		std::ifstream infile;
		oss << static_cast<std::string>(g_working_dir);
		oss << "/instance/problem/combination/quadratic_assignment/" << m_file_name << ".qap";
		infile.open(oss.str().c_str());
		if (!infile) {
			throw MyExcept("read Quadratic Assignment data error");
		}
		char *savePtr;
		while (getline(infile, Line)) {
			if (!(Keyword = strtok_s((char *)Line.c_str(), Delimiters, &savePtr)))
				continue;
			for (i = 0; i < strlen(Keyword); i++)
				Keyword[i] = toupper(Keyword[i]);
			if (!strcmp(Keyword, "OPT_OBJ")) {
				char *token = strtok_s(nullptr, Delimiters, &savePtr);
				m_optima->appendObj(std::vector<Real>(1, atof(token)));
				m_optima->setObjectiveGiven(true);
			}
			else if (!strcmp(Keyword, "OPT_SOLUTION")) {
				getline(infile, Line);
				std::istringstream ss(Line);
				VariableVector<int> temp(m_number_variables);
				for (i = 0; i < m_number_variables; ++i) {
					ss >> temp[i];
				}
				dynamic_cast<Optima<VariableVector<int>>&>(*m_optima).appendVar(temp);
				m_optima->setVariableGiven(true);
			}
		}
		infile.close();
		infile.clear();
	}

	//void QuadraticAssignment::initializeSolution(SolutionBase &s, Random *rnd) const {
	//	VariableVector<int> &x = dynamic_cast<Solution<VariableVector<int>>&>(s).variable();
	//	std::vector<int> temp;
	//	int i, pos, num = x.size();
	//	for (i = 0; i < num; i++)
	//		temp.push_back(int(i));
	//	for (i = 0; i < num; i++) {
	//		pos = int((num - 1 - i) * rnd->uniform.next());
	//		x[i] = temp[pos];
	//		temp[pos] = temp[num - 1 - i];
	//	}
	//	if (!isValid(s))
	//		throw MyExcept("error in @QuadraticAssignment::initializeSolution() in quadratic_assignment.cpp");
	//}

	void QuadraticAssignment::initializeSolution(SolutionBase &s, Random *rnd) const {
		VariableVector<int> &x = dynamic_cast<Solution<VariableVector<int>>&>(s).variable();
		std::vector<int> temp;
		int i, pos, num = x.size();
		for (i = 0; i < num; i++)
			temp.push_back(int(i));
		for (i = 0; i < num; i++) {
			pos = int((num - 1 - i) * rnd->uniform.next());
			x[i] = temp[pos];
			temp[pos] = temp[num - 1 - i];
		}
		if (!isValid(s))
			throw MyExcept("error in @QuadraticAssignment::initializeSolution() in quadratic_assignment.cpp");
	}


	bool QuadraticAssignment::same(const SolutionBase &s1, const SolutionBase &s2) const {
		const VariableVector<int> &x1 = dynamic_cast<const Solution<VariableVector<int>> &>(s1).variable();
		const VariableVector<int> &x2 = dynamic_cast<const Solution<VariableVector<int>> &>(s2).variable();
		for (int i = 0; i < m_number_variables; i++)
			if (x1[i] != x2[i])
				return false;
		return true;
	}

	Real QuadraticAssignment::variableDistance(const SolutionBase &s1, const SolutionBase &s2) const {
		const VariableVector<int> &x1 = dynamic_cast<const Solution<VariableVector<int>> &>(s1).variable();
		const VariableVector<int> &x2 = dynamic_cast<const Solution<VariableVector<int>> &>(s2).variable();
		Real dis = 0;
		for (int i = 0; i < m_number_variables; i++)
			if (x1[i] != x2[i])
				dis++;
		return dis;
	}

	Real QuadraticAssignment::variableDistance(const VariableBase &x1, const VariableBase &x2) const {
		const auto &x1_ = dynamic_cast<const VariableVector<int>&>(x1);
		const auto &x2_ = dynamic_cast<const VariableVector<int>&>(x2);
		Real dis = 0;
		for (int i = 0; i < m_number_variables; i++)
			if (x1_[i] != x2_[i])
				dis++;
		return dis;
	}

	void QuadraticAssignment::updateCandidates(const SolutionBase &sol, std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		if (candidates.empty())
			candidates.emplace_back(new Solution<VariableVector<int>>(dynamic_cast<const Solution<VariableVector<int>>&>(sol)));
		else if (sol.dominate(*candidates.front(), this))
			candidates.front().reset(new Solution<VariableVector<int>>(dynamic_cast<const Solution<VariableVector<int>>&>(sol)));
	}

	size_t QuadraticAssignment::numOptimaFound(const std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		if (m_optima->isObjectiveGiven()
			&& !candidates.empty()
			&& candidates.front()->objectiveDistance(m_optima->objective(0)) < m_objective_accuracy)
			return 1;
		else
			return 0;
	}
	
	int QuadraticAssignment::updateEvaluationTag(SolutionBase &s, Algorithm *alg) {
		if (isValid(s)) return EvaluationTag::kNormalEval;
		else return EvaluationTag::kInfeasible;
	}
}