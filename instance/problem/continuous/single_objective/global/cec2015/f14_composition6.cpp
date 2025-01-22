#include "f14_composition6.h"
#include "../../global/classical/happy_cat.h"
#include "../../global/classical/griewank_rosenbrock.h"
#include "../../global/classical/schwefel.h"
#include "../../global/classical/scaffer_F6.h"
#include "../../global/classical/elliptic.h"
#include "../../global/classical/bent_cigar.h"
#include "../../global/classical/rastrigin.h"

namespace ofec::cec2015 {
	void F14_composition6::initialize_(Environment *env) {
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

	void F14_composition6::setFunction() {
		m_num_function = 7;
		m_function.resize(m_num_function);
		m_height.resize(m_num_function);
		m_converge_severity.resize(m_num_function);
		m_f_bias.resize(m_num_function);

		BasicFunctions f(7);
		f[0] = &createFunction<HappyCat>;
		f[1] = &createFunction<GriewankRosenbrock>;
		f[2] = &createFunction<Schwefel>;
		f[3] = &createFunction<ScafferF6>;
		f[4] = &createFunction<Elliptic>;
		f[5] = &createFunction<BentCigar>;
		f[6] = &createFunction<Rastrigin>;

		for (size_t i = 0; i < m_num_function; ++i) {
			m_function[i].reset(dynamic_cast<Function *>(f[i]()));
			m_function[i]->setRandom(m_random);
			m_function[i]->setParam(m_param);
			m_function[i]->initialize();
			m_function[i]->setBias(0);
			m_function[i]->setDomain(m_domain);
			m_function[i]->setConditionNumber(1.0);
		}

		m_converge_severity[0] = 10;
		m_converge_severity[1] = 20;
		m_converge_severity[2] = 30;
		m_converge_severity[3] = 40;
		m_converge_severity[4] = 50;
		m_converge_severity[5] = 50;
		m_converge_severity[6] = 50;

		m_height[0] = 10;
		m_height[1] = 2.5;
		m_height[2] = 2.5;
		m_height[3] = 10;
		m_height[4] = 1e-6;
		m_height[5] = 1e-6;
		m_height[6] = 10;

		m_f_bias[0] = 0;
		m_f_bias[1] = 100;
		m_f_bias[2] = 200;
		m_f_bias[3] = 300;
		m_f_bias[4] = 300;
		m_f_bias[5] = 400;
		m_f_bias[6] = 400;
	}

	void F14_composition6::evaluateObjective(Real *x, std::vector<Real> &obj) {
		evaluateObjective(x, obj);
		obj[0] += 1400;
	}
}