#include "gear_train.h"

namespace ofec {


	void GearTrain::updateOptima(Environment* env)  {
		Continuous::updateOptima(env);
		Real mat_opt[4][4] = {
	{ 16,19,43,49 },{ 19,16,43,49 },{ 16,19,49,43 },{ 19,16,49,43 }
		};
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<Real> opt_obj(1);
		for (size_t i = 0; i < 4; ++i) {
			for (size_t j = 0; j < 4; ++j)
				temp_sol.variable()[j] = mat_opt[i][j];
			evaluateObjective(temp_sol.variable().vect().data(), temp_sol.objective());
			dynamic_cast<Optima<>*>(m_optima.get())->appendSolution(temp_sol);
			//m_optima->appendObj(opt_obj);
		}
		//m_optima->setVariableGiven(true);
		//std::vector<Real> opt_obj(1, 0);
		//m_optima->appendObj(opt_obj);
		//m_optima->setObjectiveGiven(true);
	}
	void GearTrain::initialize_(Environment* env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;

		resizeVariable(4);
		setDomain(12, 60);

	//	m_optima.reset(new Optima<>());

	}

	void GearTrain::evaluateObjective(Real *x, std::vector<Real> &obj) {	
		obj[0] = abs(6.931 - floor(x[2]) * floor(x[3]) / (floor(x[0]) * floor(x[1])));
	}
}