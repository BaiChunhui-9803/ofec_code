/********* Begin Register Information **********
{
    "description": "",
    "identifier": "Easom",
    "name": "Classic-Easom",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_EASOM_H
#define OFEC_EASOM_H

#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	class Easom : public Continuous {
		OFEC_CONCRETE_INSTANCE(Easom)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *vars, std::vector<Real> &objs) const override;
	};
}

#endif // ! OFEC_EASOM_H
