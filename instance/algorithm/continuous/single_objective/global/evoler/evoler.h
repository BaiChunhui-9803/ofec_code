/********* Begin Register Information **********
{
    "dependency on libraries": [ "Eigen" ],
    "identifier": "EVOLER",
    "name": "EVOLER",
    "tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/


#ifndef OFEC_EVOLER_H
#define OFEC_EVOLER_H

#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class EVOLER : virtual public Algorithm {
        OFEC_CONCRETE_INSTANCE(EVOLER)
	protected:
        size_t m_discrete_length;
        size_t m_sampling_length;

        void addInputParameters();
		void run_(Environment *env) override;

    private:
        void reconstructRepresentation(size_t M, size_t N, size_t s, Environment *env);

        std::vector<std::vector<Real>> m_grid_obj;

    public:
        const std::vector<std::vector<Real>>& gridObj() const { return m_grid_obj; }
	};
}

#endif // ! OFEC_EVOLER_H
