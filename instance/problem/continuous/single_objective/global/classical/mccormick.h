/********* Begin Register Information **********
{
    "description": "",
    "identifier": "McCormick",
    "name": "Classic-McCormick",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_MCCORMICK_H
#define OFEC_MCCORMICK_H

#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	class McCormick : public Continuous {
		OFEC_CONCRETE_INSTANCE(McCormick)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *vars, std::vector<Real> &objs) const override;
	};
}

#endif // ! OFEC_MCCORMICK_H
