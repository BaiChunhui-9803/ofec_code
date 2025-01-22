/********* Begin Register Information **********
{
	"name": "Classic_Michalewicz",
	"identifier": "Michalewicz",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_MICHALEWICZ_H
#define OFEC_MICHALEWICZ_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_mmop.h"
// modified by DYY, there is only one global optima d = 2
// https://www.sfu.ca/~ssurjano/michal.html

namespace ofec {	
	class Michalewicz : public Continuous, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(Michalewicz)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;

		int m_m;
	};
}
#endif // !OFEC_MICHALEWICZ_H