/********* Begin Register Information **********
{
	"name": "Classic_uneven_de_maxima",
	"identifier": "UnevenDeMaxima",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

//Beasley's F4 function

#ifndef OFEC_UNEVEN_DE_MAXIMA_H
#define OFEC_UNEVEN_DE_MAXIMA_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	class UnevenDeMaxima : public Continuous, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(UnevenDeMaxima)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real>& obj) const override;
	};	
}
#endif // !OFEC_UNEVEN_DE_MAXIMA_H