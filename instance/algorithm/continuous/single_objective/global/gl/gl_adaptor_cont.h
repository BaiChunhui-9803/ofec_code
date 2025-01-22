#ifndef OFEC_GL_ADAPTOR_CONT_H
#define OFEC_GL_ADAPTOR_CONT_H

#include "../../../../template/framework/gl/gl_adaptor.h"
#include "../../../../../../core/problem/solution.h"
#include <algorithm>

namespace ofec {
	class AdaptorContGL : public AdaptorGL<Solution<>> {
	public:
		AdaptorContGL(Real alpha, size_t num_dim, size_t size_pop);
		void updateProbability(Environment *env,
			PopGL<Solution<>> &pop,
			const std::vector<Real> &weight,
			const std::vector<int> *index = nullptr) override;
		void createSolution(Environment *env, Random *rnd,
			PopGL<Solution<>> &pop,
			std::vector<Solution<>> &offspring) override;
		int updateSolution(Environment *env,
			PopGL<Solution<>> &pop,
			std::vector<Solution<>> &offspring, int &num_improve) override;
		void updateStep(Environment *env, PopGL<Solution<>> &pop);
	protected:
		void accumlateProbability();
		void localSearch(Environment *env, Random *rnd, size_t i,
			PopGL<Solution<>> &pop,
			std::vector<Solution<>> &offspring);
		void globalSearch(Environment *env, Random *rnd, size_t i,
			PopGL<Solution<>> &pop,
			std::vector<Solution<>> &offspring);
	protected:
		int m_num;
		struct Infor {
			Real val;
			std::vector<int> idx;
		};
		std::vector<std::vector<Infor> > m_proba;
		std::vector<std::vector<Real> >m_acc_proba;
		struct Limit {
			std::pair<Real, Real>boundary;
			Real step, range;
			std::vector<Real> as;
		};
		std::vector<Limit> m_limit/*, m_initialRange*/;
		std::vector<std::vector<int>> m_pos;
	};
}

#endif // !OFEC_GL_ADAPTOR_CONT_H

