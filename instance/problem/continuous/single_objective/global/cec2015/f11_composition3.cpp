#include "f11_composition3.h"
#include "../classical/hg_bat.h"
#include "../classical/rastrigin.h"
#include "../classical/schwefel.h"
#include "../classical/weierstrass.h"
#include "../classical/elliptic.h"
#include "../../../../../../core/problem/solution.h"

namespace ofec::cec2015 {
	void F11_composition3::initialize_(Environment *env) {
		Continuous::initialize_();
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;

		auto& v = *m_param;
		resizeVariable(v.get<int>("number of variables"));
		setDomain(-100., 100.);
		setFunction();

		loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");
		loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");
		for (auto& i : m_function) {
			i->setGlobalOpt(i->translation().data());
		}
		// Set optimal solution
		auto new_optima = new Optima<>();
		new_optima->appendSolution(m_function[0]->optima()->solution(0));	
		m_optima.reset(new_optima);
	}

	void F11_composition3::setFunction() {
		m_num_function = 5;
		m_function.resize(m_num_function);
		m_height.resize(m_num_function);
		m_converge_severity.resize(m_num_function);
		m_f_bias.resize(m_num_function);

		BasicFunctions f(5);
		f[0] = &createFunction<HGBat>;
		f[1] = &createFunction<Rastrigin>;
		f[2] = &createFunction<Schwefel>;
		f[3] = &createFunction<Weierstrass>;
		f[4] = &createFunction<Elliptic>;

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
		m_converge_severity[2] = 10;
		m_converge_severity[3] = 20;
		m_converge_severity[4] = 20;

		m_height[0] = 10;
		m_height[1] = 10;
		m_height[2] = 2.5;
		m_height[3] = 25;
		m_height[4] = 1e-6;

		Real temp = 0;
		for (auto &i : m_f_bias) {
			i = temp;
			temp += 100;
		}
	}

	void F11_composition3::evaluateObjective(Real* x, std::vector<Real>& obj) {
		Composition::evaluateObjective(x, obj);
		obj[0] += 1100;
	}
}