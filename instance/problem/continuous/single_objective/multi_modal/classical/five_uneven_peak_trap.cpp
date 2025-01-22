#include "five_uneven_peak_trap.h"

namespace ofec {
	void FiveUnevenPeakTrap::addInputParameters() {
		m_input_parameters.add("objective accuracy", new EnumeratedReal(
			m_objective_accuracy, { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 }, 1e-4));
	}

	void FiveUnevenPeakTrap::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		resizeVariable(1);
		setDomain(0, 30);
	}

	void FiveUnevenPeakTrap::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<Real> opt_x = { 0,30 };
		for (Real x : opt_x) {
			temp_sol.variable()[0] = x;
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 0.01;
	}

	void FiveUnevenPeakTrap::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s = -1.0;
		if (x[0] >= 0 && x[0] < 2.5) {
			s = 80 * (2.5 - x[0]);
		}
		else if (x[0] >= 2.5 && x[0] < 5.0) {
			s = 64 * (x[0] - 2.5);
		}
		else if (x[0] >= 5.0 && x[0] < 7.5) {
			s = 64 * (7.5 - x[0]);
		}
		else if (x[0] >= 7.5 && x[0] < 12.5) {
			s = 28 * (x[0] - 7.5);
		}
		else if (x[0] >= 12.5 && x[0] < 17.5) {
			s = 28 * (17.5 - x[0]);
		}
		else if (x[0] >= 17.5 && x[0] < 22.5) {
			s = 32 * (x[0] - 17.5);
		}
		else if (x[0] >= 22.5 && x[0] < 27.5) {
			s = 32 * (27.5 - x[0]);
		}
		else if (x[0] >= 27.5 && x[0] <= 30) {
			s = 80 * (x[0] - 27.5);
		}

		obj[0] = s;
	}
}