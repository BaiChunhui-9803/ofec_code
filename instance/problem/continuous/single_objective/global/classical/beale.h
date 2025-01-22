/********* Begin Register Information **********
{
    "description": "",
    "identifier": "Beale",
    "name": "Classic-Beale",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_BEALE_H
#define OFEC_BEALE_H

#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	class Beale : public Continuous {
		OFEC_CONCRETE_INSTANCE(Beale)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *vars, std::vector<Real> &objs) const override;
	};
}

#endif // ! OFEC_BEALE_H
