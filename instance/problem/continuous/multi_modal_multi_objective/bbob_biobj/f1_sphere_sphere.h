/********* Begin Register Information **********
{
	"name": "BIOBJ_F01",
	"identifier": "BiobjF1",
	"problem tags": ["ConOP","GOP","MOP"]
}
*********** End Register Information **********/

/******************************************************************************************
*  Paper:2016 COCO: The Bi-objective Black Box Optimization Benchmarking (bbob-biobj)Test Suite
*******************************************************************************************/


#ifndef OFEC_F1_SPHERE_SPHERE_H
#define OFEC_F1_SPHERE_SPHERE_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"

namespace ofec {
	class BiobjF1 : public Continuous, public MetricsMOP {
	public:
		void initialize_() override;
		void evaluateObjective(Real* x, std::vector<Real>& obj);
		//void evaluateObjectiveVector(std::vector<Real>& x, std::vector<Real>& obj)
		void compute_Xopt();
		void compute_Fopt();

	};
}
#endif // 
