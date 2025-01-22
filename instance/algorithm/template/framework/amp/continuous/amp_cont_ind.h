#ifndef OFEC_AMP_Solution_CONTINUOUS_H
#define OFEC_AMP_Solution_CONTINUOUS_H


#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../core/problem/solution.h"

namespace ofec {
	template<class IndividualType>
	class IndContAMP : public IndividualType {
	protected:
		// whether the Solution fall into explored area
		bool m_explored = false;
		std::vector<Real> m_real_obj;

	public:
		IndContAMP(size_t size_obj, size_t size_con, size_t size_var) : IndividualType(size_obj, size_con, size_var) {}
		IndContAMP(const Solution<> &sol) : IndividualType(sol), m_real_obj(sol.objectiveSize()) {}
		void brwonianMove(Problem *pro, Random *rnd, double radius) {
			for (size_t j = 0; j < this->variable().size(); j++) {
                this->variable()[j] += rnd->normal.nextNonStd(0, radius);
			}
			CAST_CONOP(pro)->validateSolution(*this, Validation::kSetToBound, rnd);
		}
		void cauchyMove(Problem *pro, Random *rnd, double radius=-1) {
			auto& domain = CAST_CONOP(pro)->domain();
			for (size_t i = 0; i < this->variable().size(); i++) {
				if (radius < 0) {
                    this->variable()[i] += rnd->cauchy.nextNonStd(0, (domain[i].limit.second - domain[i].limit.first) / 2);
				}
				else {
                    this->variable()[i] += rnd->cauchy.nextNonStd(0, radius);
				}
			}
			CAST_CONOP(pro)->validateSolution(*this, Validation::kSetToBound, rnd);
		}
	};

}

#endif // ! OFEC_AMP_Solution_CONTINUOUS_H
