#include "f15_composition7.h"
#include "../../global/classical/rastrigin.h"
#include "../../global/classical/weierstrass.h"
#include "../../global/classical/happy_cat.h"
#include "../../global/classical/schwefel.h"
#include "../../global/classical/rosenbrock.h"
#include "../../global/classical/hg_bat.h"
#include "../../global/classical/katsuura.h"
#include "../../global/classical/scaffer_F6.h"
#include "../../global/classical/griewank_rosenbrock.h"
#include "../../global/classical/ackley.h"

namespace ofec::cec2015 {
	void F15_composition7::initialize_(Environment *env) {
		Continuous::initialize_();
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;

		auto &v = *m_param;
		resizeVariable(v.get<int>("number of variables"));
		setDomain(-100., 100.);
		setFunction();

		loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");
		loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");
		for (auto &i : m_function) {
			i->setGlobalOpt(i->translation().data());
		}
		// Set optimal solution
		auto new_optima = new Optima<>();
		new_optima->appendSolution(m_function[0]->optima()->solution(0));
		m_optima.reset(new_optima);
	}

	void F15_composition7::setFunction() {
		BasicFunctions f(10);
		f[0] = &createFunction<Rastrigin>;
		f[1] = &createFunction<Weierstrass>;
		f[2] = &createFunction<HappyCat>;
		f[3] = &createFunction<Schwefel>;
		f[4] = &createFunction<Rosenbrock>;
		f[5] = &createFunction<HGBat>;
		f[6] = &createFunction<katsuura>;
		f[7] = &createFunction<ScafferF6>;
		f[8] = &createFunction<GriewankRosenbrock>;
		f[9] = &createFunction<Ackley>;

		for (size_t i = 0; i < m_num_function; ++i) {
			m_function[i].reset(dynamic_cast<Function*>(f[i]()));
			m_function[i]->setRandom(m_random);
			m_function[i]->setParam(m_param);
			m_function[i]->initialize();
			m_function[i]->setBias(0);
			m_function[i]->setDomain(m_domain);
			m_function[i]->setConditionNumber(1.0);
		}

		m_converge_severity[0] = 10;
		m_converge_severity[1] = 10;
		m_converge_severity[2] = 20;
		m_converge_severity[3] = 20;
		m_converge_severity[4] = 30;
		m_converge_severity[5] = 30;
		m_converge_severity[6] = 40;
		m_converge_severity[7] = 40;
		m_converge_severity[8] = 50;
		m_converge_severity[9] = 50;


		m_height[0] = 0.1;
		m_height[1] = 2.5e-1;
		m_height[2] = 0.1;
		m_height[3] = 2.5e-2;
		m_height[4] = 1e-3;
		m_height[5] = 0.1;
		m_height[6] = 1e-5;
		m_height[7] = 10;
		m_height[8] = 2.5e-2;
		m_height[9] = 1e-3;

		m_f_bias[0] = 0;
		m_f_bias[1] = 100;
		m_f_bias[2] = 100;
		m_f_bias[3] = 200;
		m_f_bias[4] = 200;
		m_f_bias[5] = 300;
		m_f_bias[6] = 300;
		m_f_bias[7] = 400;
		m_f_bias[8] = 400;
		m_f_bias[9] = 500;
	}

	void F15_composition7::evaluateObjective(Real *x, std::vector<Real> &obj) {
		evaluateObjective(x, obj);
		obj[0] += 1500;
	}
}