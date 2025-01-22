#include "cma_es_pop.h"
#include "../../../../../../utility/functional.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "api/cmaes_interface.h"

namespace ofec {
	PopCMA_ES::PopCMA_ES(size_t size_pop, Environment *env) :
		Population(size_pop, env, env->problem()->numberVariables()), m_evo() {}

	void PopCMA_ES::resize(size_t size, Environment *env) {
		Population::resize(size, env, env->problem()->numberVariables());
	}

	void PopCMA_ES::initialize(Environment *env, Random *rnd) {
		m_ar_funvals = cmaes_init(&m_evo, env->problem()->numberVariables(), NULL, NULL, m_individuals.size(), NULL, rnd);
	}

	void PopCMA_ES::reproduce(Environment *env, Random *rnd) {
		m_pop = cmaes_SamplePopulation(&m_evo, rnd);
		for (size_t i = 0; i < m_individuals.size(); i++) {
			restoreVar(m_pop[i], m_individuals[i]->variable().vector(), env);
		}

		/* resample the solution i until it satisfies given	box constraints (variable boundaries) */
		// dead loop! by DYY

		for (size_t i = 0; i < m_individuals.size(); ++i) {
			int t = 100;
			while (CAST_CONOP(env->problem())->boundaryViolated(*m_individuals[i])&&t--) {
				cmaes_ReSampleSingle(&m_evo, i, rnd);
				restoreVar(m_pop[i], m_individuals[i]->variable().vector(), env);
			}
			// by DYY
			if (CAST_CONOP(env->problem())->boundaryViolated(*m_individuals[i])) {
				CAST_CONOP(env->problem())->validateSolution(*m_individuals[i], Validation::kSetToBound, rnd);
			}
			
		}
	}

	int PopCMA_ES::evaluate(Environment *env) {
		int tag = kNormalEval;
		m_ar_funvals = cmaes_GetArFunvals(&m_evo);
		for (size_t i = 0; i < m_individuals.size(); i++) {
			tag = m_individuals[i]->evaluate(env);
			if (tag != kNormalEval)
				break;
			m_ar_funvals[i] = mapObj(*m_individuals[i], env);
		}
		return tag;
	}

	int PopCMA_ES::evolve(Environment *env, Random *rnd) {
		cmaes_UpdateDistribution(&m_evo, m_ar_funvals);
		reproduce(env, rnd);
		auto tag = evaluate(env);
		m_iteration++;
		return tag;
	}

	PopCMA_ES::~PopCMA_ES() {
		cmaes_exit(&m_evo); /* release memory */
	}

	void PopCMA_ES::normVar(const std::vector<Real> &var, double *x, Environment *env) {
		auto &domain = CAST_CONOP(env->problem())->domain();
		for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
			x[j] = mapReal<double>(var[j], domain[j].limit.first, domain[j].limit.second, 0, 1);
		}
	}

	void PopCMA_ES::normDis(const std::vector<Real> &dis, double *s, Environment *env) {
		auto &domain = CAST_CONOP(env->problem())->domain();
		for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
			s[j] = mapReal<double>(dis[j], 0, domain[j].limit.second - domain[j].limit.first, 0, 1);
		}
	}

	void PopCMA_ES::restoreVar(double *x, std::vector<Real> &var, Environment *env) {
		auto &domain = CAST_CONOP(env->problem())->domain();
		for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
			var[j] = mapReal<double>(x[j], 0, 1, domain[j].limit.first, domain[j].limit.second);
		}
	}

	void PopCMA_ES::restoreDis(double *s, std::vector<Real> &dis, Environment *env) {
		auto &domain = CAST_CONOP(env->problem())->domain();
		for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
			dis[j] = mapReal<double>(s[j], 0, 1, 0, domain[j].limit.second - domain[j].limit.first);
		}
	}

	double PopCMA_ES::mapObj(const Solution<> &s, Environment *env) {
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			return s.objective()[0];
		}
		else {
			return -s.objective()[0];
		}
	}

	void PopCMA_ES::initializeByNonStd(Environment *env, Random *rnd, const std::vector<Real> &xstart, const std::vector<Real> &stddev) {
		std::vector<double> xstart_(xstart.size()), stddev_(stddev.size());
		normVar(xstart, xstart_.data(), env);
		normDis(stddev, stddev_.data(), env);
		m_ar_funvals = cmaes_init(&m_evo, env->problem()->numberVariables(), xstart_.data(), stddev_.data(), m_individuals.size(), NULL, rnd);
	}

	void PopCMA_ES::initializeBySample(Environment *env, Random *rnd, size_t num_samples) {
		size_t num_vars = env->problem()->numberVariables();
		VariableVector<> tmp_vars;
		std::list<std::vector<Real>> samples;
		for (size_t i = 0; i < num_samples; ++i) {
			env->problem()->initializeVariables(tmp_vars, rnd);
			samples.push_back(tmp_vars.vector());
		}
		std::vector<Real> xstart(num_vars, 0), stddev(num_vars, 0);
		for (size_t j = 0; j < num_vars; j++) {
			for (auto &sample : samples) {
				xstart[j] += sample[j];
			}
			xstart[j] /= num_samples;
			for (auto &sample : samples) {
				stddev[j] += pow(sample[j] - xstart[j], 2);
			}
			stddev[j] = sqrt(stddev[j] / (num_samples - 1));
		}
		initializeByNonStd(env, rnd, xstart, stddev);
	}

	void PopCMA_ES::initCMAES(Environment *env) {
		size_t size_pop = m_individuals.size();
		m_pop = cmaes_GetPop(&m_evo);
		for (size_t i = 0; i < size_pop; i++) {
			normVar(m_individuals[i]->variable().vector(), m_pop[i], env);
		}
		m_ar_funvals = cmaes_GetArFunvals(&m_evo);
		for (size_t i = 0; i < size_pop; i++) {
			m_ar_funvals[i] = mapObj(*m_individuals[i], env);
		}
		cmaes_InitDistribution(&m_evo, m_ar_funvals);
	}

	void PopCMA_ES::resizePopCMAES(Environment *env) {
		size_t size_pop = m_individuals.size();
		cmaes_ResizeLambda(&m_evo, size_pop);
	}

	double PopCMA_ES::conditionNumber() {
		return cmaes_Get(&m_evo, "axisratio");
	}

	bool PopCMA_ES::equalFitness() {
		return cmaes_CheckEqualFitness(&m_evo, m_ar_funvals);
	}
}
