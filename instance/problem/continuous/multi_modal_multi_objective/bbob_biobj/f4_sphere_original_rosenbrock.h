/********* Begin Register Information **********
{
	"name": "BIOBJ_F04",
	"identifier": "BiobjF4",
	"problem tags": ["MOP", "ConOP"]
}
*********** End Register Information **********/

/******************************************************************************************
*  Paper:2016 COCO: The Bi-objective Black Box Optimization Benchmarking (bbob-biobj)Test Suite
*******************************************************************************************/


#ifndef OFEC_F4_SPHERE_ORIGINAL_ROSENBROCK_H
#define OFEC_F4_SPHERE_ORIGINAL_ROSENBROCK_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"

namespace ofec {
	class BiobjF4 : public Continuous, public MetricsMOP {
	public:
		int sample_num_dimension = 70;

		void initialize_() override;
		void evaluateObjective(Real* x, std::vector<Real>& obj);
		void compute_Xopt();
		void compute_Fopt();
	};
}
#endif // 
