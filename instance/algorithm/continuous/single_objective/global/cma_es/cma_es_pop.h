#ifndef OFEC_CMAES_POP_H
#define OFEC_CMAES_POP_H

#include "../../../../../../core/algorithm/population.h"
#include "../../../../../../core/problem/solution.h"
#include "api/cmaes.h"

namespace ofec {
	class PopCMA_ES : public Population<Solution<>> {
	public:
		PopCMA_ES() = default;
		PopCMA_ES(size_t size_pop, Environment *env);
		void resize(size_t size, Environment *env);
		void initialize(Environment *env, Random *rnd) override;
		int evaluate(Environment *env) override;
		void reproduce(Environment *env, Random *rnd);
		int evolve(Environment *env, Random *rnd) override;
		~PopCMA_ES();
		void initializeByNonStd(Environment *env, Random *rnd, const std::vector<Real> &xstart, const std::vector<Real> &stddev);
		void initializeBySample(Environment *env, Random *rnd, size_t num_samples = 1000);
		void initCMAES(Environment *env);
		void resizePopCMAES(Environment *env);
		double conditionNumber();
		bool equalFitness();

	protected:
		cmaes_t m_evo;
		double *m_ar_funvals = nullptr, *const *m_pop = nullptr, *m_xfinal = nullptr;

	protected:
		void normVar(const std::vector<Real> &var, double *x, Environment *env);
		void normDis(const std::vector<Real> &dis, double *s, Environment *env);
		void restoreVar(double *x, std::vector<Real> &var, Environment *env);
		void restoreDis(double *s, std::vector<Real> &dis, Environment *env);
		double mapObj(const Solution<> &s, Environment *env);
	};
}

#endif // !OFEC_CMAES_POP_H