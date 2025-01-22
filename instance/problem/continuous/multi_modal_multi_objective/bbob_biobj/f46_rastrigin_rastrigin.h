/********* Begin Register Information **********
{
	"name": "BIOBJ_F46",
	"identifier": "BiobjF46",
	"problem tags": ["ConOP","GOP","MOP"]
}
*********** End Register Information **********/

/******************************************************************************************
*  Paper:2016 COCO: The Bi-objective Black Box Optimization Benchmarking (bbob-biobj)Test Suite
*******************************************************************************************/


#ifndef OFEC_F46_RASTRIGIN_RASTRIGIN_H
#define OFEC_F46_RASTRIGIN_RASTRIGIN_H

#include "../../../../../core/problem/continuous/function.h"
#include "../../../../../utility/linear_algebra/matrix.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"

namespace ofec {
	class BiobjF46 : public Continuous, public MetricsMOP {
	public:
		int sample_num_dimension = 120;
		Matrix m_rot, m_rot2, m_linearTF;
		Real m_beta = 0.2;
		Real m_condition_number;
		std::vector<Real> m_norRand;

		void initialize_() override;
		void evaluateObjective(Real* x, std::vector<Real>& obj);
		void compute_Xopt();
		void compute_Fopt();

		void irregularize(Real* x);
		void asyemmetricalize(Real* x, Real belta);
		void reshape(Matrix& B, std::vector<Real>& vector, size_t m, size_t n);
		void computeRotation(Matrix& rot, size_t Dim);
		bool loadRotation(Real base_);


	};
}
#endif // 
