/********* Begin Register Information **********
{
	"name": "RS-CMSA",
	"identifier": "RS_CMSA",
	"tags": [ "continuous", "single-objective" ],
	"dependency on libraries": [ "Eigen" ]
}
*********** End Register Information **********/

/*
  Ahrari, A., Deb, K., & Preuss, M. (2016). 
  Multimodal Optimization by Covariance Matrix Self-Adaptation Evolution Strategy with Repelling Subpopulations. 
  Evolutionary computation, doi:10.1162/EVCO_a_00182
*/

#ifndef OFEC_RS_CMSA_H
#define OFEC_RS_CMSA_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "api/struct.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../core/algorithm/population.h"
#include "../../../../../../core/problem/solution.h"

namespace ofec {
	class RS_CMSA : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(RS_CMSA)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

	private:
		size_t m_remained_eval;		// Evaluation budget
		size_t m_targetNsubP;		// Default number of subpopulations
		rs_cmsa::Opt m_opt;		// Algorithm parameters
		MultiPopulation<Population<Solution<>>> m_subpops;
	};
}

#endif // !OFEC_RS_CMSA_H
