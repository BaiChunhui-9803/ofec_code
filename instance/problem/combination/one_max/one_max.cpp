#include "one_max.h"
#include "../../../../core/problem/solution.h"
#include <map>

namespace ofec {

	void OneMax::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1e6, 100));
	}

	void OneMax::initialize_(Environment* env) {
		ProblemVariableVector::initialize_(env);
		resizeObjective(1);
		resizeConstraint(0);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
	}

	void OneMax::updateOptima(Environment* env) {
		m_optima.reset(new Optima<VariableType>());
		SolutionType curOpt(m_number_objectives, m_number_constraints, m_number_variables);
		std::fill(curOpt.variable().vect().begin(), curOpt.variable().vect().end(), 1);
		curOpt.objective()[0] = m_number_variables;
		dynamic_cast<Optima<VariableType>&>(*m_optima).appendSolution(curOpt);
	}

	void OneMax::evaluate(const VariableBase& vars, std::vector<Real>& objs, std::vector<Real>& cons) {
		auto &x = dynamic_cast<const VariableType&>(vars);
		for (int i = 0; i < m_number_objectives; i++)
			objs[i] = 0;
		for (int n = 0; n < m_number_objectives; n++)
			for (size_t i = 0; i < m_number_variables; i++)
				objs[n] += x[i];
	}


	bool OneMax::isValid(const SolutionBase &s) {
		if (!m_if_valid_check) return true;
		const auto &x = dynamic_cast<const Solution<VariableVector<int>> &>(s).variable();
		for (int i = 0; i < m_number_variables; i++) {
			if (x[i] != 0 && x[i] != 1)
				return false;
		}
		return true;
	}

	void OneMax::initializeVariables(VariableBase& s, Random *rnd) const
	{
		auto& x = dynamic_cast<VariableType&>(s);
		for (size_t i = 0; i < m_number_variables; i++)
			if (rnd->uniform.next() < 0.5)
				x[i] = 0;
			else x[i] = 1;
	}

	
	bool OneMax::same(const VariableBase&s1, const VariableBase&s2) const {
		const auto &x1 = dynamic_cast<const Solution<VariableVector<int>> &>(s1).variable();
		const auto &x2 = dynamic_cast<const Solution<VariableVector<int>> &>(s2).variable();
		for (int i = 0; i < m_number_variables; i++)
			if (x1[i] != x2[i])
				return false;
		return true;
	}

	Real OneMax::variableDistance(const VariableBase&x1, const VariableBase&x2) const {
		const auto& x1_ = dynamic_cast<const VariableType&>(x1);
		const auto& x2_ = dynamic_cast<const VariableType&>(x2);
		Real dis = 0;
		for (int i = 0; i < m_number_variables; i++)
			if (x1_[i] != x2_[i])
				dis++;
		return dis;
	}

	void OneMax::initSolutionNearBy(const SolutionBase& centerSol, SolutionBase& sol, int neighbork, ofec::Random* rnd) {
		const auto& centerX = dynamic_cast<const Solution<VariableVector<int>> &>(centerSol).variable();
		auto& solX = dynamic_cast<Solution<VariableVector<int>> &>(sol).variable();
		solX = centerX;
		for (int idx(0); idx < neighbork; ++idx) {
			solX[rnd->uniform.nextNonStd<int>(0, m_number_variables)] = 1 - solX[rnd->uniform.nextNonStd<int>(0, m_number_variables)];
		}
		
	}

	void OneMax::initSolVecNearBy(const SolutionBase& centerSol, SolutionBase& sol, std::vector<bool>& vec) {
		const auto& centerX = dynamic_cast<const Solution<VariableVector<int>> &>(centerSol).variable();
		auto& solX = dynamic_cast<Solution<VariableVector<int>> &>(sol).variable();
		solX = centerX;
		for (int idx(0); idx < m_number_variables; ++idx) {
			if (vec[idx]) {
				solX[idx] = 1 - solX[idx];
			}
		}
	}


	void OneMax::filterSameSols(std::vector<const SolutionBase*>& sols, 
		std::vector<bool>& uniqueFlag) {

		uniqueFlag.resize(sols.size());
		std::fill(uniqueFlag.begin(), uniqueFlag.end(), true);

		ofec::Random rnd(0.5);
		std::vector<std::array<unsigned, 2>> randNum(m_number_variables);
		std::array<unsigned, 2> min_max = { std::numeric_limits<unsigned>::max()/2, std::numeric_limits<unsigned>::max() };
		for (auto& it : randNum) {
			for (auto& it2 : it) {
				it2 = rnd.uniform.nextNonStd<unsigned>(min_max.front(), min_max.back());
			}
		}

		std::map<unsigned, std::set<std::vector<int>>> id2sols;

		for (int idSol(0);idSol<sols.size();++idSol) {
			auto& cursol = sols[idSol];
			auto& solX = dynamic_cast<const Solution<VariableVector<int>> &>(*cursol).variable().vect();
			unsigned curRandNum(0);
			for (int idx(0); idx < m_number_variables; ++idx) {
				curRandNum ^= randNum[idx][solX[idx]];
			}

			if (id2sols.find(curRandNum) == id2sols.end()) {
				id2sols[curRandNum].insert(solX);
			}
			else {
				auto& mapSols = id2sols[curRandNum];
				int originSize = mapSols.size();
				mapSols.insert(solX);
				if (originSize == mapSols.size()) {
					uniqueFlag[idSol] = false;
				}
			}
		}
	}
}


