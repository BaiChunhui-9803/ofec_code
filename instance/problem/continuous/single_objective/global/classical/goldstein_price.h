/********* Begin Register Information **********
{
    "description": "",
    "identifier": "GoldsteinPrice",
    "name": "Classic-Goldstein-Price",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_GOLDSTEIN_PRICE_H
#define OFEC_GOLDSTEIN_PRICE_H

#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	class GoldsteinPrice : public Continuous {
		OFEC_CONCRETE_INSTANCE(GoldsteinPrice)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *vars, std::vector<Real> &objs) const override;
	};
}

#endif // ! OFEC_GOLDSTEIN_PRICE_H
