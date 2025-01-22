#ifndef OFEC_RS_CMSA_ES_ANALYZE_CONVERGED_SOL_H
#define OFEC_RS_CMSA_ES_ANALYZE_CONVERGED_SOL_H

#include "../../../../../../../utility/linear_algebra/eigen.h"
#include "struct.h"
#include "../../../../../../../core/environment/environment.h"

namespace ofec {
	class Problem;
	class Algorithm;

	namespace rs_cmsa {
		using namespace Eigen;

		void analyzeConvergedSol(const Ref<const MatrixXd> &Y, Environment *env, MatrixXd &archiveXf, VectorXd &all_N_rep,
			VectorXd &NormTabDis, const Opt &opt, size_t remained_eval, size_t &used_eval);
	}
}

#endif //!OFEC_RS_CMSA_ES_ANALYZE_CONVERGED_SOL_H